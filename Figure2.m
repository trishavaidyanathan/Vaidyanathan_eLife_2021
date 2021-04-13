%% %%%%% FIGURE 2 %%%%%%

%% Data sets used in this figure
% Endogenous Ca WT
pathWT = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\Endogenous Ca2+\WT'; %change this to your path on your computer
% Endogenous Ca KO
pathKO = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\Endogenous Ca2+\KO'; %change this to your path on your computer

%Load endogenous Ca2+ (WT): 
load(strcat(pathWT,'\ephysData'))
load(strcat(pathWT,'\aquaData'))

%Load endogenous Ca2+ (KO): 
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
         
%% Figure 2A 

%%% Using WT data: 

sf = 1000; %sampling frequency
for d = 5 %For example, picked d = 5 181016 m111 004, 2200:3100 s
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
    ha(1) = subplot(2,1,1);
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
    ha(2) = subplot(2,1,2);
        histogram(allData(d).peakTimes,1:length(allData(d).numEvents))
    
    linkaxes(ha,'x')
    set(gcf,'renderer','painters')
end
xlim([2200 3100])

%% Figure 2B
%%% Uses the spect, T, F output from Fig2A ^

d = 5; %For example, picked d = 5 181016 m111 004, 2200:3100 s
% the time windows chosen for each example trace:
locoWake = 2372:2382; 
statWake = 2352:2362;
sleep = 2800:2810;
locoWakeF = round(2372/allData(d).frameperiod):round(2382/allData(d).frameperiod); 
statWakeF = round(2352/allData(d).frameperiod):round(2362/allData(d).frameperiod);
sleepF = round(2800/allData(d).frameperiod):round(2810/allData(d).frameperiod);


%%% CHANGE THIS TO THE BEHAVIORAL PERIOD YOU WANT TO PLOT:
rangeF = sleepF; range = sleep; name = 'sleep';
%%%

[binValues, binEdges] = histcounts(allData(d).peakTimesF,1:allData(d).numFrames);
binEdges = binEdges(binEdges >=rangeF(1) & binEdges <= rangeF(end));
binValues = binValues(binEdges);

gap = min([1/max(binValues),.25]);
spikeTimes = [];
for i = 1:length(binEdges)
    xs = (([1:binValues(i)]-(binValues(i)/2))*gap)+(binEdges(i)+.5);    
    spikeTimes = [spikeTimes xs];
end

figure;
    ha(1) = subplot(2,1,1);
        plot(T, mean(spect,2),'--r','linewidth',1.5)    
        ylim([.015 .04])
        title(name)
        xlim([range(1) range(end)])
    ha(2) = subplot(2,1,2);
        rasterplot(spikeTimes-rangeF(1),1,length(rangeF),ha(2))
        ylim([-1 2])
        hold on;
        vline(1:length(rangeF),'-k')
        xlim([1 length(rangeF)]); xlabel('frames')

%% Figure 2C 

%%% WT data should be loaded as "allDataWT" and "allResWT" AND
%%% KO data should be loaded as "allDataIP3R2KO" and "allResIP3R2KO"
%%% Uses function "CaETSpect.m", make sure it's in your path!

behav = {'locomotory wake', 'sleep','stationary wake'};
window = [5 5];
datastruct = allDataWT;
resstruct = allResWT;
preDeltaWT = {}; postDeltaWT = {}; spectWT = {};
for b = 1:3 %behavioral periods: wake, sleep, stat wake
    events = {};
    for i = 1:length(datastruct)
        if b == 1
        events{i} = [datastruct(i).movingEvtIdx{:}];
        elseif b == 2
        events{i} = [datastruct(i).sleepEvtIdx{:}];
        elseif b == 3
        events{i} = [datastruct(i).notmovingEvtIdx{:}];
        end
    end
    name = strcat('WT',{' '},behav{b});
    [pre, post, Spect, ~,T,F] = CaETSpect(datastruct, resstruct, events, window, lfp, name);
    preDeltaWT{b} = pre;
    postDeltaWT{b} = post;
    spectWT{b} = Spect;
end

behav = {'locomotory wake', 'sleep','stationary wake'};
window = [5 5];
datastruct = allDataIP3R2KO;
resstruct = allResIP3R2KO;
preDeltaKO = {}; postDeltaKO = {}; spectKO = {}; 
for b = 1:3 %behavioral periods: wake, sleep, stat wake
    events = {};
    for i = 1:length(datastruct)
        if b == 1
        events{i} = [datastruct(i).movingEvtIdx{:}];
        elseif b == 2
        events{i} = [datastruct(i).sleepEvtIdx{:}];
        elseif b == 3
        events{i} = [datastruct(i).notmovingEvtIdx{:}];
        end
    end
    disp(strcat('n = ',num2str(sum(cellfun(@length,events)))))
    name = strcat('IP3R2KO',{' '},behav{b});
    [pre, post, Spect, ~,T,F] = CaETSpect(datastruct, resstruct, events,window, lfp, name);
    preDeltaKO{b} = pre;
    postDeltaKO{b} = post;
    spectKO{b} = Spect;
end


% Fig 2C LEFT (WT)
%Plot shaded error bar
colors = {'-b','-r','-c'}; clear a;
figure; 
for b = 2:3
    deltaTraceMean = squeeze(median(spectWT{b},3)); 
    deltaTraceMean = mean(deltaTraceMean,2);
    
    deltaTrace = squeeze(mean(spectWT{b},2));
    deltaTraceSEM = std(deltaTrace,[],2) ./ size(deltaTrace,2);
    
    s = shadedErrorBar(T,deltaTraceMean,deltaTraceSEM,'lineprops',colors{b},'transparent',1); 
    a(b) = s.mainLine;
    hold on; 
end
xlim([1 9]); ylim([.94 1.04]);
v = vline(window(1),'-k'); v.LineWidth = 2;
h = hline(1,'--k'); h.LineWidth = 2;
xlabel('sec'); ylabel(strcat('normalized power:',{' '},'0.5-4 Hz'));
title('WT: ETA Astrocyte Event Onsets')
legend(a(2:3),behav(2:3),'location','northwest')
set(gcf,'renderer','painters')

% Fig 2C RIGHT (KO)
%shaded error bar for KO
colors = {'-b','-r','-c'}; clear a;
figure; 
for b = 2:3
    deltaTraceMean = squeeze(median(spectKO{b},3));  
    deltaTraceMean = mean(deltaTraceMean,2);
    
    deltaTrace = squeeze(mean(spectKO{b},2));
    deltaTraceSEM = std(deltaTrace,[],2) ./ size(deltaTrace,2);
    
    s = shadedErrorBar(T,deltaTraceMean,deltaTraceSEM,'lineprops',colors{b},'transparent',1); 
    a(b) = s.mainLine;
    hold on; 
end
xlim([1 9]); ylim([.94 1.04]);
v = vline(window(1),'-k'); v.LineWidth = 2;
h = hline(1,'--k'); h.LineWidth = 2;
xlabel('sec'); ylabel(strcat('normalized power:',{' '},'0.5-4 Hz'));
title('ETA Astrocyte Event Onsets: IP3R2KO')
legend(a(2:3),behav(2:3),'location','northwest')
set(gcf,'renderer','painters')

%% Figure 2D

PostPreDiffWT = cellfun(@minus, postDeltaWT, preDeltaWT, 'UniformOutput', 0);
meansWT = cellfun(@mean, PostPreDiffWT);
SEMsWT = cellfun(@sem, PostPreDiffWT);

PostPreDiffKO = cellfun(@minus, postDeltaKO, preDeltaKO, 'UniformOutput', 0);
meansKO = cellfun(@mean, PostPreDiffKO);
SEMsKO = cellfun(@sem, PostPreDiffKO);

figure;
s1 =scatter([1 2],[meansWT(2) meansKO(2)],50,'filled');
hold on;
errorbar([1 2],[meansWT(2) meansKO(2)],[SEMsWT(2) SEMsKO(2)],'k','linestyle','none');
hold on;
s2 =scatter([1.2 2.2],[meansWT(3) meansKO(3)],50,'filled');
hold on;
errorbar([1.2 2.2],[meansWT(3) meansKO(3)],[SEMsWT(3) SEMsKO(3)],'k','linestyle','none');
[p12_WT] = (ranksum(PostPreDiffWT{2}, PostPreDiffWT{3}));
[p12_KO] = (ranksum(PostPreDiffKO{2}, PostPreDiffKO{3}));
[~,p11] = (ttest2(PostPreDiffWT{2}, PostPreDiffKO{2}));
[~,p22] = (ttest2(PostPreDiffWT{3}, PostPreDiffKO{3}));
hold on;
textmax = double(max(meansWT+SEMsWT));
text(1,textmax,strcat('WT',':',{' '},num2str(p12_WT)))
text(2,textmax,strcat('KO',':',{' '},num2str(p12_KO)))
text(1.5,textmax,strcat('sleep',':',{' '},num2str(p11)))
text(1.5,textmax-.01,strcat('statWake',':',{' '},num2str(p22)))
xticks([1 2]); xticklabels({'WT','IP3R2KO'});
ylabel('mean delta after onset - mean delta before onset')
title('mean delta 2 sec before vs mean delta 2 sec after')
xlim([.5 2.5]); 
hold on;
hline(0,'--k')
set(gca,'fontsize',15)
legend([s1 s2],{'sleep','stationary wake'},'location','southwest')
