%% %%%%% FIGURE 1 FIGURE SUPPLEMENT 1 %%%%%%

%% Data sets used in this figure
% Endogenous Ca WT
pathWT = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\Endogenous Ca2+\WT'; %change this to your path on your computer
% PCA Endogenous Ca 
pathPCA = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\PCA\Endogenous Ca2+'; %change this to your path on your computer

%PCA Output 
load(strcat(pathPCA,'\PCA_Output_Endogenous'));

%For endogenous Ca2+ (WT) (Panel A): 
load(strcat(pathWT,'\ephysData'))
load(strcat(pathWT,'\aquaData'))

%% Pre-processing FOR ENDOGENOUS Ca2+ DATASET
    %% Artifact removal
        %%%% Not sure what files I did what to exactly

        %filter out super-low frequencies for drifting baseline
        for d = []
            order = 10000;
            cutoffFreq = .3;

            highpassFilt = designfilt('highpassfir', 'FilterOrder', order, 'CutoffFrequency', cutoffFreq, 'SampleRate', 1000); %designfilter

            % fvtool(highpassFilt) %view filter

            new = filtfilt(highpassFilt, double(allData(d).baseline_lfp)); %use filter on LFP data 
            figure; plot(new)
            allData(d).baseline_lfp = new;
        end

        %Remove artifacts by SD thresh
        for d = []; 
            toremove = find(zscore(allData(d).lfp) > 6 | zscore(allData(d).lfp) < -5.5 );
            figure; plot(allData(d).times, allData(d).lfp); hold on;
            scatter(allData(d).times(toremove),allData(d).lfp(toremove),'filled')

            allData(d).lfp(toremove) = mean(allData(d).lfp) ;
        end 
    %% Create uniqueMice cell-array
        sf=1000;

        [uniqueMice, ia, ic] = unique({allData(:).mouseID}); %unqueMice contains all mouseID strings
        %ia contains the index number of allData corresponding to first instance of the mouse ID
        %ic contains the index number of uniqueMice for each entry of allData 
        mice = {}; %each entry corresponds to uniqueMice, contains the indeces of allData for each mouse
        for m = 1:length(uniqueMice) %for each mouse, find indeces corresponding to CNO, Saline, CNO-hour1, Saline-hour1, etc
            currMouse = find(ic == m);
            mice{m} = currMouse; %each entry of mice has a vector of indeces of allData for that mouse       
        end
    %% SLEEP-SCORE (Mult-files)
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

%% Fig1FigSupp1 A

legendNames = {'Wake-locomotion', 'Wake-stationary', 'SWS'};
titleNames = {'event rate (log)', 'amplitude', 'amplitude mean', 'size (log)', 'size mean (log)', 'duration', 'duration mean', 'delta power (log)'};
% titleNames = {'event rate', 'amplitude', 'amplitude mean', 'size (log)', 'size mean (log)', 'duration', 'duration mean', 'delta power (log)'};
pngNames = {'event_rate', 'amplitude', 'amplitude_mean', 'size', 'size_mean', 'duration', 'duration_mean','delta'};

eventRateAll = {log([allData(:).movingEventRate]), log([allData(:).notmovingEventRate]), log([allData(:).sleepEventRate])};
%eventRateAll = {([allData(:).movingEventRate]), ([allData(:).notmovingEventRate]), ([allData(:).sleepEventRate])};
DffMaxAll = {[allData(:).movingDffMax], [allData(:).notmovingDffMax], [allData(:).sleepDffMax]};
DffPeriodAll = {[allData(:).movingDffMax_PeriodMean], [allData(:).notmovingDffMax_PeriodMean], [allData(:).sleepDffMax_PeriodMean]};
AreaAll = {log([allData(:).movingArea]), log([allData(:).notmovingArea]), log([allData(:).sleepArea])};
AreaPeriodAll = {log([allData(:).movingArea_PeriodMean]), log([allData(:).notmovingArea_PeriodMean]), log([allData(:).sleepArea_PeriodMean])};
DurationAll = {[allData(:).movingDuration], [allData(:).notmovingDuration], [allData(:).sleepDuration]};
DurationPeriodAll = {[allData(:).movingDuration_PeriodMean], [allData(:).notmovingDuration_PeriodMean], [allData(:).sleepDuration_PeriodMean]};
deltaPeriod = {[allData(:).movingDelta_PeriodMean], [allData(:).notmovingDelta_PeriodMean], [allData(:).sleepDelta_PeriodMean]};


all = {eventRateAll DffMaxAll DffPeriodAll AreaAll AreaPeriodAll DurationAll DurationPeriodAll deltaPeriod};

clear toPlot
for i = 1:length(titleNames)
    toPlot(i).titleName = titleNames{i};
    toPlot(i).pngName = pngNames{i};
    
    curr = all{i};
    toPlot(i).sleepData = (curr{3});    
    toPlot(i).sleepData(isinf(toPlot(i).sleepData)) = 0;
    toPlot(i).statWakeData = (curr{2});
    toPlot(i).statWakeData(isinf(toPlot(i).statWakeData)) = 0;
    toPlot(i).locoWakeData = (curr{1});
    toPlot(i).locoWakeData(isinf(toPlot(i).locoWakeData)) = 0;      
end

for i = [5 7 3]
    curr = toPlot(i);
    figure;
        violinplot(curr)%,[],'ShowData',false)
        p12 = (ranksum(curr.sleepData,curr.statWakeData));
        p13 = (ranksum(curr.sleepData,curr.locoWakeData));
        p23 = (ranksum(curr.statWakeData,curr.locoWakeData));
        hold on;
        textmax = max([curr.locoWakeData curr.statWakeData curr.sleepData]);
        textmin = min([curr.locoWakeData curr.statWakeData curr.sleepData]);   
        text(.5,textmax,strcat('1 and 2',':',{' '},num2str(p12)))
        text(1.5,textmax+1,strcat('1 and 3',':',{' '},num2str(p13)))
        text(2,textmax,strcat('2 and 3',':',{' '},num2str(p23)))
        ylim([textmin-1 textmax+2])
        ylabel(titleNames{i});
    title(titleNames{i})
    set(gca,'fontsize',15)    
    set(gcf,'renderer','painters')    
end

%% Fig1FigSupp1 B

names = {'area (um^2)',...
    'perimeter (um)',...
    'circularity',...
    'amplitude',...
    'amplitude with exclusion',...
    'rise time (s)',...
    'fall time (s)',...
'duration 50% width (s)', ...
'duration 10% width (s)', ...
'propagation growth (anterior)',...
'propagation growth (posterior)',...
'propagation growth (left)',...
'propagation growth (right)',...
'propagation shrink (anterior)',...
'propagation shrink (posterior)',...
'propagation shrink (left)',...
'propagation shrink (right)',...
'temporal overlap',...
'spatial overlap',...
'spatial overlap size corrected'}; %From featureNames 

order = [1:3, 19:20, 10:17, 6:9, 18, 4:5]; 

figure;
subplot(4,1,1:3)
imagesc(weights(:,order)')
yticks([1:size(weights,2)]); yticklabels(names(order));
xticks([1:size(weights,1)]); xticklabels({'PC1','PC2','PC3','PC4','PC5'})
colorbar;
set(gca,'fontsize',15)
subplot(4,1,4)
percentExplained = explainedVarianceRatio*100;
bar(percentExplained); colorbar; ylabel('% explained');
ylim([0 60])
text(1:length(explainedVarianceRatio),percentExplained,strsplit(num2str(percentExplained)),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom'); 
set(gca,'fontsize',15)

%% Remaining figures in Fig1FigSupp1 were done in Python
