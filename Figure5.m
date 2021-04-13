%% %%%%% FIGURE 5 %%%%%%

%% Data sets used in this figure
% Endogenous Ca WT
pathWT = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\Endogenous Ca2+\WT'; %change this to your path on your computer
% Endogenous Ca KO
pathKO = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\Endogenous Ca2+\KO'; %change this to your path on your computer

%For endogenous Ca2+ (WT)
load(strcat(pathWT,'\ephysData'))
load(strcat(pathWT,'\aquaData'))

%For endogenous Ca2+ (KO)
load(strcat(pathKO,'\ephysData'))
load(strcat(pathKO,'\aquaData'))

%% Pre-processing (Assumes data is in "allData" and "allRes") 
    %% Optional: Artifact removal

        %filter out super-low frequencies for drifting baseline
        for d = [] %fill in depending on what files need to be filetered
            order = 10000;
            cutoffFreq = .3;

            highpassFilt = designfilt('highpassfir', 'FilterOrder', order, 'CutoffFrequency', cutoffFreq, 'SampleRate', 1000); %designfilter

            % fvtool(highpassFilt) %view filter

            new = filtfilt(highpassFilt, double(allData(d).baseline_lfp)); %use filter on LFP data 
            figure; plot(new)
            allData(d).baseline_lfp = new;
        end

        %Remove artifacts by SD thresh
        for d = []
            toremove = find(zscore(allData(d).lfp) > 6 | zscore(allData(d).lfp) < -5.5 );
            figure; plot(allData(d).times, allData(d).lfp); hold on;
            scatter(allData(d).times(toremove),allData(d).lfp(toremove),'filled')

            allData(d).lfp(toremove) = mean(allData(d).lfp);      
        end
    %% Create uniqueMice cell-array
        sf=1000; %sampling frequency

        [uniqueMice, ia, ic] = unique({allData(:).mouseID}); %unqueMice contains all mouseID strings
        %ia contains the index number of allData corresponding to first instance of the mouse ID
        %ic contains the index number of uniqueMice for each entry of allData 
        mice = {}; %each entry corresponds to uniqueMice, contains the indeces of allData for each mouse
        for m = 1:length(uniqueMice) %for each mouse, find indeces corresponding to CNO, Saline, CNO-hour1, Saline-hour1, etc
            currMouse = find(ic == m);
            mice{m} = currMouse; %each entry of mice has a vector of indeces of allData for that mouse       
        end
    %% SLEEP-SCORE (Mult-files)
        % Requries Chronux! 
        sf = 1000; %samples per sec
        tic
        for m = 1:length(uniqueMice)
            files = mice{m};
            data = {}; names = {};
            k=1;
            for f = files'
                temp = [allData(f).times, allData(f).lfp, allData(f).emg, round(allData(f).opto)];
                data{k} = temp;
                names{k} = allData(f).name;
                k = k+1;
            end

            [INT, IDX, f_FFT, FFT, t_FFT, v_FFT, swTimes] = ...
                sleepDetection_zscore_multFiles(data, sf, names, 0);

            k = 1;
            for f = files'
                allData(f).INT = (INT{k});
                allData(f).IDX = single(IDX{k});
                allData(f).f_FFT = single(f_FFT{k});
                allData(f).FFT = single(FFT{k});
                allData(f).t_FFT = single(t_FFT{k});
                allData(f).v_FFT = single(v_FFT{k});
                allData(f).swTimes = single(swTimes{k});
                k=k+1;
            end

            disp(strcat('done with mouse number',{' '},uniqueMice{m}));
        end
        toc
    %% Find spectral info
        freqList = [.5 4 ;
            6 12 ;
            15 30 ;
            40 60 ;
            65 100];

        % Each row of freqList is a frequency range corresponding to freqNames

        freqNames = {'delta', 'theta', 'beta', 'low gamma', 'high gamma'};
        % freqNames = {'delta', 'theta', 'beta', 'gamma'};

        %Iterate through freq bands and find average power over time using LFP
        for d = 1:length(allData)
            specSelect = []; %each row is power over time, corresponds to freqList
            specSelectNorm = []; %contians normalized power over time (divide by sum of all other powers)
            for i = 1:size(freqList,1)
                first = freqList(i,1);
                last = freqList(i,2);

                [~,f_first] = min(abs(allData(d).f_FFT-first));
                [~,f_last] = min(abs(allData(d).f_FFT-last));

                temp = log(allData(d).FFT(:,f_first:f_last));
                specSelect(i,:) = mean(temp,2);       
            end
            allData(d).specSelect = specSelect;
        end
    %% Calculate locomotion
        win = 1000; %window in ms to look at opto
        movementBufferWin = 10; %in seconds: how many sec before start of stationary period must there be no locomotion

        for d = 1:length(allData)

            winstart = [0:win:(length(allData(d).times)-win)]+1; %index corresponding to the start of each window

            allData(d).locomotion = []; %vector of 0 (not moving) and 1 (moving) for each window (1 sec)
            for i = 1:length(winstart) 
                m = winstart(i);
                optowin = allData(d).opto(m:m+(win-1));
                optowinDif = diff(optowin); %will be all 0 if every value in optowin is the same
                if isempty(find(optowinDif,1)) %if optowinDif are all 0's
                    allData(d).locomotion(i) = 0; %0 means not moving
                else
                    allData(d).locomotion(i) = 1; %1 = moving
                end
            end
            allData(d).whenMoving = []; %each row is a moving period, 1st column is start of movement period, C2 is end of movement period in sec (like INT but for movement)    
            allData(d).whenNotMoving = [];
            for i = 1:size(allData(d).INT{1,1},1)
                start = allData(d).INT{1,1}(i,1); %when wake period starts in sec
                finish = allData(d).INT{1,1}(i,2); %finish of wake period in sec
                start = round(start);
                finish = round(finish);

                temp = find(allData(d).locomotion(start:finish)) + (start-1); %finds indices of locomotion where locomotion is 1 (ie the sec where there is movement)
                if ~isempty(temp)
                    move = find(diff(temp)>1); %for multiple mvmt periods in temp, this finds the ends of each one (besides the last one)
                    e = [temp(move), temp(end)]; %marks all the end points of the movement period
                    s = [temp(1) temp(move + 1)]; %marks all the start points of the movement period

                    allData(d).whenMoving = [allData(d).whenMoving ; s' e'];
                end

                temp = find(allData(d).locomotion(start:finish) == 0) + (start-1); %finds idx where locomotion is 0 (ie the sec where there is movement)
                if ~isempty(temp)
                    Nomove = find(diff(temp)>1); %when there are >1 period, marks the ends of each stationary period (besidse the last one)
                    e = [temp(Nomove), temp(end)]; %marks all the end points of the movement period
                    s = [temp(1) temp(Nomove + 1)]; %marks all the start points of the movement period

            %Add criteria that this stationary wake period must be at least min duration
                    lengths = e-s;
                    toKeep = lengths > 15;

                    s = s(toKeep);
                    e = e(toKeep);

            %Add condition here where non-moving periods cannot immediately preceed a moving period - must have a buffer (movementBufferWin)                    
                    for j = 1:length(s)
                        if allData(d).locomotion(s(j)-1) == 1 %if start of stationary period immediately follows movement
                            s(j) = s(j) + movementBufferWin;
                        end
                    end

                    allData(d).whenNotMoving = [allData(d).whenNotMoving ; s' e'];
                end
            end
        end
    %% Calculate event rate
        for d = 1:length(allData)
            allData(d).frameperiod = allRes(d).opts.frameRate;
            allData(d).numFrames = allRes(d).opts.sz(3);
            allData(d).peakTimesF = allRes(d).fts.loc.t0; %in frames (onsets)
            allData(d).peakTimes = allData(d).peakTimesF*allData(d).frameperiod; %in sec
            [allData(d).numEvents, allData(d).timebins] = hist(allData(d).peakTimes, 1:(allData(d).numFrames * allData(d).frameperiod)); %
            allData(d).eventRate = allData(d).numEvents/(allData(d).timebins(2)-allData(d).timebins(1)); %events per second
        end

        for d = 1:length(allData)
            [allData(d).onsetSec, allData(d).EvtsSortByOnset] = sort(allData(d).peakTimes); %EvtsSortByOnset is list of event IDs in order of when they happen
        end
    %% Analyze based on behavioral state: sleep    

        minSWS = 15; %in sec: min duration of sleep period to be used

        for d = 1:length(allData)
            sleepNumEvents = [];
            sleepEventRate = [];
            sleepDffMax = [];
            sleepDffMax_PeriodMean = []; %each entry is the mean dff for that sleep period
            sleepArea = []; %area in pixels
            sleepArea_PeriodMean = []; %each entry is the mean size (in pixels) for that sleep period
            sleepDuration = []; %duration of events in frames (peak 10% to peak 10%)
            sleepDuration_PeriodMean = []; %each entry is mean duration (in frames) for that sleep period

            sleepDelta_PeriodMean = []; %each entry is mean delta power for that sleep period

            allData(d).sleepEvtIdx = {};

            for i = 1:size(allData(d).INT{1,2},1)
                start = allData(d).INT{1,2}(i,1); %when sleep period starts in sec
                finish = allData(d).INT{1,2}(i,2); %when sleep period ends in sec

                if finish-start > minSWS

                    [~, startIdx] = min(abs(allData(d).timebins - start)); %find idx of timebins corresponding to start of sleep period
                    [~, finishIdx] = min(abs(allData(d).timebins - finish)); %find idx of timebins corresponding to end of sleep period

                    sleepNumEvents = [sleepNumEvents sum(allData(d).numEvents(startIdx:finishIdx))];
                    sleepEventRate = [sleepEventRate mean(allData(d).eventRate(startIdx:finishIdx))];

                    relEvents = allData(d).EvtsSortByOnset(find(allData(d).onsetSec >= start & allData(d).onsetSec <= finish)); %event IDs that occur in this period
                    allData(d).sleepEvtIdx{i} = relEvents; 

                    sleepDffMax = [sleepDffMax allRes(d).fts.curve.dffMax(relEvents)];
                    sleepDffMax_PeriodMean = [sleepDffMax_PeriodMean mean(allRes(d).fts.curve.dffMax(relEvents))];
                    sleepArea = [sleepArea allRes(d).fts.basic.area(relEvents)];
                    sleepArea_PeriodMean = [sleepArea_PeriodMean mean(allRes(d).fts.basic.area(relEvents))];
                    sleepDuration = [sleepDuration allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents)];
                    sleepDuration_PeriodMean = [sleepDuration_PeriodMean mean(allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents))];


                    %For Delta/Theta
                    [~, startIdx] = min(abs(start-allData(d).t_FFT)); %find start of sleep period using output of lfp analysis (t)
                    [~, finishIdx] = min(abs(finish-allData(d).t_FFT)); %find end of sleep period using output of lfp analysis

                    sleepDelta_PeriodMean = [sleepDelta_PeriodMean, mean(allData(d).specSelect(1,startIdx:finishIdx))];
                end
            end

            allData(d).sleepNumEvents = sleepNumEvents;
            allData(d).sleepEventRate = sleepEventRate;
            allData(d).sleepDffMax = sleepDffMax;
            allData(d).sleepDffMax_PeriodMean = sleepDffMax_PeriodMean; %each entry is the mean dff for that sleep period
            allData(d).sleepArea = sleepArea; %area in pixels
            allData(d).sleepArea_PeriodMean = sleepArea_PeriodMean; %each entry is the mean size (in pixels) for that sleep period
            allData(d).sleepDuration = sleepDuration; %duration of events in frames (peak 10% to peak 10%)
            allData(d).sleepDuration_PeriodMean = sleepDuration_PeriodMean; %each entry is mean duration (in frames) for that sleep period

            allData(d).sleepDelta_PeriodMean = sleepDelta_PeriodMean; %each entry is mean delta power for that sleep period
        end
    %% Analyze based on behavioral state: wake   

        for d = 1:length(allData)
            wakeNumEvents = [];
            wakeEventRate = [];
            wakeDffMax = [];
            wakeDffMax_PeriodMean = []; %each entry is the mean dff for that wake period
            wakeArea = []; %area in pixels
            wakeArea_PeriodMean = []; %each entry is the mean size (in pixels) for that wake period
            wakeDuration = []; %duration of events in frames (peak 10% to peak 10%)
            wakeDuration_PeriodMean = []; %each entry is mean duration (in frames) for that wake period

            wakeDelta_PeriodMean = []; %each entry is mean delta power for that sleep period

            allData(d).wakeEvtIdx = {};

            for i = 1:size(allData(d).INT{1,1},1)
                start = allData(d).INT{1,1}(i,1); %when sleep period starts in sec
                finish = allData(d).INT{1,1}(i,2); %when sleep period ends in sec

                [~, startIdx] = min(abs(allData(d).timebins - start)); %find idx of timebins corresponding to start of sleep period
                [~, finishIdx] = min(abs(allData(d).timebins - finish)); %find idx of timebins corresponding to end of sleep period

                wakeNumEvents = [wakeNumEvents sum(allData(d).numEvents(startIdx:finishIdx))];
                wakeEventRate = [wakeEventRate mean(allData(d).eventRate(startIdx:finishIdx))];

                relEvents = allData(d).EvtsSortByOnset(allData(d).onsetSec >= start & allData(d).onsetSec <= finish); %event IDs that occur in this period
                allData(d).wakeEvtIdx{i} = relEvents;

                wakeDffMax = [wakeDffMax allRes(d).fts.curve.dffMax(relEvents)];
                wakeDffMax_PeriodMean = [wakeDffMax_PeriodMean mean(allRes(d).fts.curve.dffMax(relEvents))];
                wakeArea = [wakeArea allRes(d).fts.basic.area(relEvents)];
                wakeArea_PeriodMean = [wakeArea_PeriodMean mean(allRes(d).fts.basic.area(relEvents))];
                wakeDuration = [wakeDuration allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents)];
                wakeDuration_PeriodMean = [wakeDuration_PeriodMean mean(allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents))];

                %For Delta/Theta
                [~, startIdx] = min(abs(start-allData(d).t_FFT)); %find start of sleep period using output of lfp analysis (t)
                [~, finishIdx] = min(abs(finish-allData(d).t_FFT)); %find end of sleep period using output of lfp analysis

                wakeDelta_PeriodMean = [wakeDelta_PeriodMean, mean(allData(d).specSelect(1,startIdx:finishIdx))];
            end

            allData(d).wakeNumEvents = wakeNumEvents;
            allData(d).wakeEventRate = wakeEventRate;
            allData(d).wakeDffMax = wakeDffMax;
            allData(d).wakeDffMax_PeriodMean = wakeDffMax_PeriodMean; %each entry is the mean dff for that wake period
            allData(d).wakeArea = wakeArea; %area in pixels
            allData(d).wakeArea_PeriodMean = wakeArea_PeriodMean; %each entry is the mean size (in pixels) for that wake period
            allData(d).wakeDuration = wakeDuration; %duration of events in frames (peak 10% to peak 10%)
            allData(d).wakeDuration_PeriodMean = wakeDuration_PeriodMean; %each entry is mean duration (in frames) for that wake period

            allData(d).wakeDelta_PeriodMean = wakeDelta_PeriodMean; %each entry is mean delta power for that sleep period
        end
    %% Analyze based on behavioral state: wake-moving 

        for d = 1:length(allData)
            movingNumEvents = [];
            movingEventRate = [];
            movingDffMax = [];
            movingDffMax_PeriodMean = []; %each entry is the mean dff for that wake period
            movingArea = []; %area in pixels
            movingArea_PeriodMean = []; %each entry is the mean size (in pixels) for that wake period
            movingDuration = []; %duration of events in frames (peak 10% to peak 10%)
            movingDuration_PeriodMean = []; %each entry is mean duration (in frames) for that wake period

            movingDelta_PeriodMean = []; %each entry is mean delta power for that sleep period

            allData(d).movingEvtIdx = {}; 

            for i = 1:length(allData(d).whenMoving)
                start = allData(d).whenMoving(i,1); %when moving period starts in sec
                finish = allData(d).whenMoving(i,2); %when moving period ends in sec

                [~, startIdx] = min(abs(allData(d).timebins - start)); %find idx of timebins corresponding to start of sleep period
                [~, finishIdx] = min(abs(allData(d).timebins - finish)); %find idx of timebins corresponding to end of sleep period

                movingNumEvents = [movingNumEvents sum(allData(d).numEvents(startIdx:finishIdx))];
                movingEventRate = [movingEventRate mean(allData(d).eventRate(startIdx:finishIdx))];

                relEvents = allData(d).EvtsSortByOnset(find(allData(d).onsetSec >= start & allData(d).onsetSec <= finish)); %event IDs that occur in this period
                allData(d).movingEvtIdx{i} = relEvents;

                movingDffMax = [movingDffMax allRes(d).fts.curve.dffMax(relEvents)];
                movingDffMax_PeriodMean = [movingDffMax_PeriodMean mean(allRes(d).fts.curve.dffMax(relEvents))];
                movingArea = [movingArea allRes(d).fts.basic.area(relEvents)];
                movingArea_PeriodMean = [movingArea_PeriodMean mean(allRes(d).fts.basic.area(relEvents))];
                movingDuration = [movingDuration allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents)];
                movingDuration_PeriodMean = [movingDuration_PeriodMean mean(allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents))];

                %For Delta/Theta
                [~, startIdx] = min(abs(start-allData(d).t_FFT)); %find start of sleep period using output of lfp analysis (t)
                [~, finishIdx] = min(abs(finish-allData(d).t_FFT)); %find end of sleep period using output of lfp analysis

                movingDelta_PeriodMean = [movingDelta_PeriodMean, mean(allData(d).specSelect(1,startIdx:finishIdx))];
            end

            allData(d).movingNumEvents = movingNumEvents;
            allData(d).movingEventRate = movingEventRate;
            allData(d).movingDffMax = movingDffMax;
            allData(d).movingDffMax_PeriodMean = movingDffMax_PeriodMean; %each entry is the mean dff for that wake period
            allData(d).movingArea = movingArea; %area in pixels
            allData(d).movingArea_PeriodMean = movingArea_PeriodMean; %each entry is the mean size (in pixels) for that wake period
            allData(d).movingDuration = movingDuration; %duration of events in frames (peak 10% to peak 10%)
            allData(d).movingDuration_PeriodMean = movingDuration_PeriodMean; %each entry is mean duration (in frames) for that wake period

            allData(d).movingDelta_PeriodMean = movingDelta_PeriodMean; %each entry is mean delta power for that sleep period
        end
    %% Analyze based on behavioral state: wake-stationary

        minStationary = 15; %in sec: min duration of a stationary period
        for d = 1:length(allData) %for each data file
            %initialize vectors:
            notmovingNumEvents = [];
            notmovingEventRate = []; %each entry is the event rate for each stationary period in this file
            notmovingDffMax = [];
            notmovingDffMax_PeriodMean = []; %each entry is the mean dff for that wake period
            notmovingArea = []; %area in pixels
            notmovingArea_PeriodMean = []; %each entry is the mean size (in pixels) for that wake period
            notmovingDuration = []; %duration of events in frames (peak 10% to peak 10%)
            notmovingDuration_PeriodMean = []; %each entry is mean duration (in frames) for that wake period
            notmovingDelta_PeriodMean = []; %each entry is mean delta power for that sleep period

            allData(d).notmovingEvtIdx = {};

            for i = 1:length(allData(d).whenNotMoving) %for each stationary period within this data file
                start = allData(d).whenNotMoving(i,1); %when stationary period starts in sec
                finish = allData(d).whenNotMoving(i,2); %when stationary period ends in sec

                if finish-start > minStationary

                    [~, startIdx] = min(abs(allData(d).timebins - start)); %find idx of timebins corresponding to start of sleep period
                    [~, finishIdx] = min(abs(allData(d).timebins - finish)); %find idx of timebins corresponding to end of sleep period

                    notmovingNumEvents = [notmovingNumEvents sum(allData(d).numEvents(startIdx:finishIdx))];
                    notmovingEventRate = [notmovingEventRate mean(allData(d).eventRate(startIdx:finishIdx))];

                    relEvents = allData(d).EvtsSortByOnset(find(allData(d).onsetSec >= start & allData(d).onsetSec <= finish)); %event IDs that occur in this period

                    allData(d).notmovingEvtIdx{i} = relEvents;

                    notmovingDffMax = [notmovingDffMax allRes(d).fts.curve.dffMax(relEvents)];
                    notmovingDffMax_PeriodMean = [notmovingDffMax_PeriodMean mean(allRes(d).fts.curve.dffMax(relEvents))];
                    notmovingArea = [notmovingArea allRes(d).fts.basic.area(relEvents)];
                    notmovingArea_PeriodMean = [notmovingArea_PeriodMean mean(allRes(d).fts.basic.area(relEvents))];
                    notmovingDuration = [notmovingDuration allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents)];
                    notmovingDuration_PeriodMean = [notmovingDuration_PeriodMean mean(allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents))];

                    %For Delta/Theta
                    [~, startIdx] = min(abs(start-allData(d).t_FFT)); %find start of sleep period using output of lfp analysis (t)
                    [~, finishIdx] = min(abs(finish-allData(d).t_FFT)); %find end of sleep period using output of lfp analysis

                    notmovingDelta_PeriodMean = [notmovingDelta_PeriodMean, mean(allData(d).specSelect(1,startIdx:finishIdx))];
                end
            end

            allData(d).notmovingNumEvents = notmovingNumEvents;
            allData(d).notmovingEventRate = notmovingEventRate;
            allData(d).notmovingDffMax = notmovingDffMax;
            allData(d).notmovingDffMax_PeriodMean = notmovingDffMax_PeriodMean; %each entry is the mean dff for that wake period
            allData(d).notmovingArea = notmovingArea; %area in pixels
            allData(d).notmovingArea_PeriodMean = notmovingArea_PeriodMean; %each entry is the mean size (in pixels) for that wake period
            allData(d).notmovingDuration = notmovingDuration; %duration of events in frames (peak 10% to peak 10%)
            allData(d).notmovingDuration_PeriodMean = notmovingDuration_PeriodMean; %each entry is mean duration (in frames) for that wake period

            allData(d).notmovingDelta_PeriodMean = notmovingDelta_PeriodMean; %each entry is mean delta power for that sleep period
        end

%% Figure 5A  
%%% LOAD WT DATA 

 %%%% Cumulative event count leading up to SWS-W transition
    win = 30; % how many sec of win to look at for each transition
    minDur = win*2; %only include sleep/wake periods that are at least twice the window of interest length

    evtCount_SWStransW = [];
    for d = 1:length(allData)

        %first find SWStransW times
        wLengths = allData(d).INT{1}(:,2)-allData(d).INT{1}(:,1);
        sLengths = allData(d).INT{2}(:,2)-allData(d).INT{2}(:,1);

        wakeStarts = allData(d).INT{1}(:,1);
        wakeStarts = wakeStarts(wLengths>minDur);
        
        sleepEnds = allData(d).INT{2}(:,2);
        sleepStarts = allData(d).INT{2}(:,1); %also find sleep starts to mark preceeding sleep period
        sleepEnds = sleepEnds(sLengths>minDur);
        sleepStarts = sleepStarts(sLengths>minDur);
        
        SWStransW = sleepEnds(ismember(sleepEnds,wakeStarts));
        sleepStarts = sleepStarts(ismember(sleepEnds,wakeStarts));
        sleeps = [sleepStarts SWStransW]; %marks start and end of sleep periods corresponding to relevant transitions

        curr_numEvts = [];
         for i = 1:length(SWStransW)
             currTrans = round(SWStransW(i));
             currSstart = round(sleeps(i,1));
             if (currTrans-win)>0
                 numEvts = allData(d).numEvents(currTrans-win:currTrans);
                 totalSevents = sum(allData(d).numEvents(currSstart:currTrans)); %use to find percent of events in that period in the win sec before transition
                 curr_numEvts = [curr_numEvts ; (numEvts/totalSevents)*100];    
             end
         end
         if size(curr_numEvts,1)>1 %for cases where there were no valid transitions or only 1
             evtCount_SWStransW= [evtCount_SWStransW ; curr_numEvts];
         end
    end

  %%%% Cumulative event count leading up to W-SWS transition
    win = 30; % how many sec of win to look at for each transition
    minDur = win*2; %only include sleep/wake periods that are at least twice the window of interest length

    evtCount_WtransSWS = [];
    for d = 1:length(allData)

        %first find WtransSWS times
        wLengths = allData(d).INT{1}(:,2)-allData(d).INT{1}(:,1);
        sLengths = allData(d).INT{2}(:,2)-allData(d).INT{2}(:,1);

        wakeEnds = allData(d).INT{1}(:,2);
        wakeStarts = allData(d).INT{1}(:,1);
        wakeEnds = wakeEnds(wLengths>minDur);
        wakeStarts = wakeStarts(wLengths>minDur);
        
        sleepStarts = allData(d).INT{2}(:,1);
        sleepStarts = sleepStarts(sLengths>minDur);
        
        WtransSWS = wakeEnds(ismember(wakeEnds,sleepStarts));
        wakeStarts = wakeStarts(ismember(wakeEnds,sleepStarts));
        wakes = [wakeStarts WtransSWS]; %marks start and end of wake periods correspodning to WtransSWS

        curr_numEvts = []; 
         for i = 1:length(WtransSWS)
             currTrans = round(WtransSWS(i));
             currWstart = round(wakes(i,1));
             if (currTrans-win)>0
                 numEvts = allData(d).numEvents(currTrans-win:currTrans);
                 totalWevents = sum(allData(d).numEvents(currWstart:currTrans)); %use to find percent of events in that period in the win sec before transition             
                 curr_numEvts = [curr_numEvts ; (numEvts / totalWevents)*100];
             end
         end
         if size(curr_numEvts,1)>1 %for cases where there were no valid transitions or only 1
             evtCount_WtransSWS = [evtCount_WtransSWS ; curr_numEvts];
         end
    end     
    
  %%%% Plot
    figure;
        xs = [fliplr(0:win)]*-1;
        ys1 = [(evtCount_SWStransW)];
        ys2 = [(evtCount_WtransSWS)];
        semCurr1 = []; for k = 1:size(ys1,2); semCurr1 = [semCurr1 nanstd(ys1(:,k))/sqrt(length(ys1(:,k)))]; end
        semCurr2 = []; for k = 1:size(ys2,2); semCurr2 = [semCurr2 nanstd(ys2(:,k))/sqrt(length(ys2(:,k)))]; end
        a = shadedErrorBar(xs,nanmean(ys1),semCurr1,'lineprops','-r','transparent',1);
        hold on;
        b = shadedErrorBar(xs,nanmean(ys2),semCurr2,'lineprops','-c','transparent',1);
        xlabel('sec before transition'); ylabel('% of total events in period');
        legend([a.mainLine b.mainLine],{'SWS -> W', 'W -> SWS'})
        title(strcat(num2str(win),[' ','sec before transition (min dur= ',' '],num2str(minDur),[' ','sec)']))
        set(gca,'fontsize',15);

%% Figure 5B     
%%% LOAD IP3R2KO DATA as "allData" and "allRes" AND PREPROCESS %%%


 %%%% Cumulative event count leading up to SWS-W transition
    win = 30; % how many sec of win to look at for each transition
    minDur = win*2; %only include sleep/wake periods that are at least twice the window of interest length

    evtCount_SWStransW = [];
    for d = 1:length(allData)

        %first find SWStransW times
        wLengths = allData(d).INT{1}(:,2)-allData(d).INT{1}(:,1);
        sLengths = allData(d).INT{2}(:,2)-allData(d).INT{2}(:,1);

        wakeStarts = allData(d).INT{1}(:,1);
        wakeStarts = wakeStarts(wLengths>minDur);
        
        sleepEnds = allData(d).INT{2}(:,2);
        sleepStarts = allData(d).INT{2}(:,1); %also find sleep starts to mark preceeding sleep period
        sleepEnds = sleepEnds(sLengths>minDur);
        sleepStarts = sleepStarts(sLengths>minDur);
        
        SWStransW = sleepEnds(ismember(sleepEnds,wakeStarts));
        sleepStarts = sleepStarts(ismember(sleepEnds,wakeStarts));
        sleeps = [sleepStarts SWStransW]; %marks start and end of sleep periods corresponding to relevant transitions

        curr_numEvts = [];
         for i = 1:length(SWStransW)
             currTrans = round(SWStransW(i));
             currSstart = round(sleeps(i,1));
             if (currTrans-win)>0
                 numEvts = allData(d).numEvents(currTrans-win:currTrans);
                 totalSevents = sum(allData(d).numEvents(currSstart:currTrans)); %use to find percent of events in that period in the win sec before transition
                 curr_numEvts = [curr_numEvts ; (numEvts/totalSevents)*100];    
             end
         end
         if size(curr_numEvts,1)>1 %for cases where there were no valid transitions or only 1
             evtCount_SWStransW= [evtCount_SWStransW ; curr_numEvts];
         end
    end

  %%%% Cumulative event count leading up to W-SWS transition
    win = 30; % how many sec of win to look at for each transition
    minDur = win*2; %only include sleep/wake periods that are at least twice the window of interest length

    evtCount_WtransSWS = [];
    for d = 1:length(allData)

        %first find WtransSWS times
        wLengths = allData(d).INT{1}(:,2)-allData(d).INT{1}(:,1);
        sLengths = allData(d).INT{2}(:,2)-allData(d).INT{2}(:,1);

        wakeEnds = allData(d).INT{1}(:,2);
        wakeStarts = allData(d).INT{1}(:,1);
        wakeEnds = wakeEnds(wLengths>minDur);
        wakeStarts = wakeStarts(wLengths>minDur);
        
        sleepStarts = allData(d).INT{2}(:,1);
        sleepStarts = sleepStarts(sLengths>minDur);
        
        WtransSWS = wakeEnds(ismember(wakeEnds,sleepStarts));
        wakeStarts = wakeStarts(ismember(wakeEnds,sleepStarts));
        wakes = [wakeStarts WtransSWS]; %marks start and end of wake periods correspodning to WtransSWS

        curr_numEvts = []; 
         for i = 1:length(WtransSWS)
             currTrans = round(WtransSWS(i));
             currWstart = round(wakes(i,1));
             if (currTrans-win)>0
                 numEvts = allData(d).numEvents(currTrans-win:currTrans);
                 totalWevents = sum(allData(d).numEvents(currWstart:currTrans)); %use to find percent of events in that period in the win sec before transition             
                 curr_numEvts = [curr_numEvts ; (numEvts / totalWevents)*100];
             end
         end
         if size(curr_numEvts,1)>1 %for cases where there were no valid transitions or only 1
             evtCount_WtransSWS = [evtCount_WtransSWS ; curr_numEvts];
         end
    end     
    
  %%%% Plot
    figure;
        xs = [fliplr(0:win)]*-1;
        ys1 = [(evtCount_SWStransW)];
        ys2 = [(evtCount_WtransSWS)];
        semCurr1 = []; for k = 1:size(ys1,2); semCurr1 = [semCurr1 nanstd(ys1(:,k))/sqrt(length(ys1(:,k)))]; end
        semCurr2 = []; for k = 1:size(ys2,2); semCurr2 = [semCurr2 nanstd(ys2(:,k))/sqrt(length(ys2(:,k)))]; end
        a = shadedErrorBar(xs,nanmean(ys1),semCurr1,'lineprops','-r','transparent',1);
        hold on;
        b = shadedErrorBar(xs,nanmean(ys2),semCurr2,'lineprops','-c','transparent',1);
        xlabel('sec before transition'); ylabel('% of total events in period');
        legend([a.mainLine b.mainLine],{'SWS -> W', 'W -> SWS'})
        title(strcat(num2str(win),[' ','sec before transition (min dur= ',' '],num2str(minDur),[' ','sec)']))
        set(gca,'fontsize',15);

%% Figure 5C

%%% Load different MAT file that has m113 006 11-01-18:
load('E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\Endogenous Ca2+\WT\schematicData.mat')

sf = 1000;
d = 1; %For example, picked d = 1 181101 m113 008, 900:1500 s
    fpass = [.1,4];
    lfp = allData(d).lfp;
    % highpass
    [b,a] = butter(2,fpass(1)/(sf/2),'high');
    lfp_delta = filtfilt(b,a,lfp);

    % lowpass
    [b,a] = butter(4,fpass(2)/(sf/2),'low');
    lfp_delta = filtfilt(b,a,lfp_delta);

    movingwin = [2,0.1]; %moving window for spectrogram
    %LFP params structure
    Params.Fs = 1000;
    Params.fpass = [.1 4];
    Params.pad =-1;
    Params.tapers = [2.5 2];
    Params.trialave = 0;

    [S,T,F] = mtspecgramc(allData(d).lfp,movingwin,Params);
    spect = mean(S,2);
    
    figure;
        plot([1:length(lfp_delta)]/sf,lfp_delta,'k');
        hold on;
        for n=1:size(allData(d).INT{2},1)
            s = allData(d).INT{2}(n,1);
            fin = allData(d).INT{2}(n,2);
            h = fill([s s fin fin], [-.5 .5 .5 -.5], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            hold on;
        end
        for n=1:size(allData(d).whenMoving,1)
            s = allData(d).whenMoving(n,1);
            fin = allData(d).whenMoving(n,2);
            h = fill([s s fin fin], [-.5 .5 .5 -.5], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            hold on;
        end
        for n=1:size(allData(d).whenNotMoving,1)
            s = allData(d).whenNotMoving(n,1);
            fin = allData(d).whenNotMoving(n,2);
            h = fill([s s fin fin], [-.5 .5 .5 -.5], 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            hold on;
        end
xlim([900 1500])

%% Figure 5D

%%% LOAD WT DATA AS "allDataWT" and "allResWT" AND
%%% LOAD IP3R2KO DATA AS "allDataIP3R2KO' and "allResIP3R2KO


minDur = 60;
n = 3; %divide each sleep/wake period into n bins

sleepNumEvts_WT = []; wakeNumEvts_WT = [];  %each row is a sleep period, each column is a bin (n bins)
for d = 1:length(allDataWT)    
    sleeps = allDataWT(d).INT{2};
    wakes = allDataWT(d).INT{1};
    sLengths = sleeps(:,2)-sleeps(:,1);
    wLengths = wakes(:,2)-wakes(:,1);
    sleeps = sleeps(sLengths>minDur,:); %select sleeps that meet min duration
    wakes = wakes(wLengths>minDur,:); %select wakes that meet min duration   
    for s = 1:size(sleeps,1) %for each sleep period
        st = sleeps(s,1);
        fn = sleeps(s,2);
        bins = linspace(st,fn,n+1); 
        binSt = bins(1:end-1); %marks the start of each bin
        binFn = bins(2:end); %marks end of each bin
        numEvts = []; 
        for b = 1:length(binSt) %for each bin
            currSt = round(binSt(b));
            currFn = round(binFn(b));
            numEvts = [numEvts sum(allDataWT(d).numEvents(currSt:currFn))];
        end
        sleepNumEvts_WT = [sleepNumEvts_WT ; numEvts];
    end
    
    for w = 1:size(wakes,1) %for each sleep period
        st = wakes(w,1);
        fn = wakes(w,2);
        bins = linspace(st,fn,n+1); 
        binSt = bins(1:end-1); %marks the start of each bin
        binFn = bins(2:end); %marks end of each bin
        numEvts = [];
        for b = 1:length(binSt) %for each bin
            currSt = round(binSt(b));
            currFn = round(binFn(b));
            numEvts = [numEvts sum(allDataWT(d).numEvents(currSt:currFn))];
        end
        wakeNumEvts_WT = [wakeNumEvts_WT ; numEvts]; 
   end        
end

sleepNumEvts_KO = []; wakeNumEvts_KO = [];  %each row is a sleep period, each column is a bin (n bins)
for d = 1:length(allDataIP3R2KO)    
    sleeps = allDataIP3R2KO(d).INT{2};
    wakes = allDataIP3R2KO(d).INT{1};
    sLengths = sleeps(:,2)-sleeps(:,1);
    wLengths = wakes(:,2)-wakes(:,1);
    sleeps = sleeps(sLengths>minDur,:); %select sleeps that meet min duration
    wakes = wakes(wLengths>minDur,:); %select wakes that meet min duration   
    for s = 1:size(sleeps,1) %for each sleep period
        st = sleeps(s,1);
        fn = sleeps(s,2);
        bins = linspace(st,fn,n+1); 
        binSt = bins(1:end-1); %marks the start of each bin
        binFn = bins(2:end); %marks end of each bin
        numEvts = []; 
        for b = 1:length(binSt) %for each bin
            currSt = round(binSt(b));
            currFn = round(binFn(b));
            numEvts = [numEvts sum(allDataIP3R2KO(d).numEvents(currSt:currFn))];
        end
        sleepNumEvts_KO = [sleepNumEvts_KO ; numEvts];
    end
    
    for w = 1:size(wakes,1) %for each sleep period
        st = wakes(w,1);
        fn = wakes(w,2);
        bins = linspace(st,fn,n+1); 
        binSt = bins(1:end-1); %marks the start of each bin
        binFn = bins(2:end); %marks end of each bin
        numEvts = []; 
        for b = 1:length(binSt) %for each bin
            currSt = round(binSt(b));
            currFn = round(binFn(b));
            numEvts = [numEvts sum(allDataIP3R2KO(d).numEvents(currSt:currFn))];
        end
        wakeNumEvts_KO = [wakeNumEvts_KO ; numEvts]; 
   end        
end

%%% Figure 5D Left

figure;
    sleepCurrWT = sleepNumEvts_WT - sleepNumEvts_WT(:,1);
    sleepSEMWT = nanstd(sleepCurrWT) / sqrt(size(sleepCurrWT,1));
    e = errorbar(1:n, nanmean(sleepCurrWT), sleepSEMWT);
    e.Marker = 'o'; e.MarkerSize = 8;
    e.Color = 'r'; e.MarkerFaceColor = 'r';
hold on
    sleepCurrKO = sleepNumEvts_KO - sleepNumEvts_KO(:,1);
    sleepSEMKO = nanstd(sleepCurrKO) / sqrt(size(sleepCurrKO,1));
    e = errorbar(1:n, nanmean(sleepCurrKO), sleepSEMKO);
    e.Marker = 'o'; e.MarkerSize = 8;
    e.Color = 'm'; e.MarkerFaceColor = 'm';
hold on
xticks(1:n); 
ylabel('change in number of events')
xlim([.5 n+.5])    
set(gca,'fontsize',15);

for i = 1:n
    top = max([nanmean(sleepCurrWT(:,i))+sleepSEMWT(i), nanmean(sleepCurrKO(:,i))+sleepSEMKO(i)])+1;
    [h,p] = ttest2(sleepCurrWT(:,i), sleepCurrKO(:,i))
    if p < 0.05
        text(i,top,'*','FontSize',24);
    else
        text(i,top,'n.s.','FontSize',14);
    end
end
[h,p] = ttest(sleepCurrWT(:,i), sleepCurrWT(:,1));   
text(2.5, 2.5, strcat('WT:',{' '},'3 vs 1:',' p =',{' '},num2str(p)));
[h,p] = ttest(sleepCurrKO(:,i), sleepCurrKO(:,1));   
text(2.5, 2, strcat('KO:',{' '},'3 vs 1:',' p =',{' '},num2str(p)));
hline(0,'--k');
legend({'sleep-WT','sleep-KO'},'location','southwest')
title('sleep')

%%% Figure 5D Right

figure;
    wakeCurrWT = wakeNumEvts_WT - wakeNumEvts_WT(:,1);
    wakeSEMWT = nanstd(wakeCurrWT) / sqrt(size(wakeCurrWT,1));
    k = errorbar(1:n, nanmean(wakeCurrWT), wakeSEMWT);
    k.Marker = 'o'; k.MarkerSize = 8;
    k.Color = 'b'; k.MarkerFaceColor = 'b';
hold on
    wakeCurrKO = wakeNumEvts_KO - wakeNumEvts_KO(:,1);
    wakeSEMKO = nanstd(wakeCurrWT) / sqrt(size(wakeCurrWT,1));
    k = errorbar(1:n, nanmean(wakeCurrKO), wakeSEMKO);
    k.Marker = 'o'; k.MarkerSize = 8;
    k.Color = 'c'; k.MarkerFaceColor = 'c';
xticks(1:n); 
ylabel('change in number of events')
xlim([.5 n+.5])    
set(gca,'fontsize',15);
for i = 1:n
    top = max([nanmean(wakeCurrWT(:,i))+wakeSEMWT(i),nanmean(wakeCurrKO(:,i))+wakeSEMKO(i)])+1;
    [h,p] = ttest2(wakeCurrWT(:,i), wakeCurrKO(:,i));
    if p < 0.05
        text(i,top,'*','FontSize',24);
    else
        text(i,top,'n.s.','FontSize',14);
    end
end
[h,p] = ttest(wakeCurrWT(:,i), wakeCurrWT(:,1));   
text(2.5, -10, strcat('WT:',{' '},'3 vs 1:',' p =',{' '},num2str(p)));
[h,p] = ttest(wakeCurrKO(:,i), wakeCurrKO(:,1));   
text(2.5, -13, strcat('KO:',{' '},'3 vs 1:',' p =',{' '},num2str(p)));
hline(0,'--k');
legend({'wake-WT','wake-KO'},'location','southwest')
title('wake')
