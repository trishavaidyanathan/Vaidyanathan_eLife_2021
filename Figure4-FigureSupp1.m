%% %%%%% FIGURE 4 FIGURE SUPPLEMENT 1 %%%%%%

%% Data sets used in this figure
% Gi DREADD WT
pathWT = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\DREADDs\Gi DREADDs\WT'; %change this to your path on your computer

%Load Gi-DREADD WT 
load(strcat(pathWT,'\ephysData'))
load(strcat(pathWT,'\aquaData'))

%% Pre-processing (Assumes data is in "allData" and "allRes") 
    %% Create uniqueMice cell-array
        sf=1000; %sampling frequency

        [uniqueMice, ia, ic] = unique({allData(:).mouseID}); %unqueMice contains all mouseID strings
        %ia contains the index number of allData corresponding to first instance of the mouse ID
        %ic contains the index number of uniqueMice for each entry of allData 
        mice = {}; %each entry corresponds to uniqueMice, contains the indeces of allData for each mouse
        CNO = {}; Sal = {}; %each entry corresponds to uniqueMice, contains the indeces of allData for all CNO/Sal files only of that mouse
        CNO_1 = {}; CNO_2 = {}; CNO_3 = {}; %each entry corresponds to uniqueMice, contains the indeces of allData for just CNO hour1/2/3 files
        Sal_1 = {}; Sal_2 = {}; Sal_3 = {}; %each entry corresponds to uniqueMice, contains the indeces of allData for just Sal hour1/2/3 files
        CNO_baseline = {};
        sal_baseline = {};
        Gq = []; Gi = []; ip3r2KO = [];
        for m = 1:length(uniqueMice) %for each mouse, find indeces corresponding to CNO, Saline, CNO-hour1, Saline-hour1, etc
            currMouse = find(ic == m);
            mice{m} = currMouse; %each entry of mice has a vector of indeces of allData for that mouse

            if ~cellfun(@isempty, (strfind({allData(currMouse).dreadd}, 'IP3R2KO')))
                ip3r2KO = [ip3r2KO m];
            elseif ~cellfun(@isempty, (strfind({allData(currMouse).dreadd}, 'Gq')))
                Gq = [Gq m];
            else
                Gi = [Gi m];
            end

            currCNObaseline = find(~cellfun(@isempty, strfind({allData(currMouse).condition}, 'CNO'))&...
                ~cellfun(@isempty, strfind({allData(currMouse).condition}, 'baseline')));
            currCNObaseline = currMouse(currCNObaseline);
            currSalbaseline = find(~cellfun(@isempty, strfind({allData(currMouse).condition}, 'Saline'))&...
                ~cellfun(@isempty, strfind({allData(currMouse).condition}, 'baseline')));
            currSalbaseline = currMouse(currSalbaseline); 
            CNO_baseline{m} = currCNObaseline;
            sal_baseline{m} = currSalbaseline;

            currCNO = find(~cellfun(@isempty, strfind({allData(currMouse).condition}, 'CNO'))&...
                cellfun(@isempty, strfind({allData(currMouse).condition}, 'baseline')));
            currCNO = currMouse(currCNO);
            currSal = find(~cellfun(@isempty, strfind({allData(currMouse).condition}, 'Saline'))&...
                cellfun(@isempty, strfind({allData(currMouse).condition}, 'baseline')));
            currSal = currMouse(currSal);
            CNO{m} = currCNO; %each entry corresponds to uniqueMice and contains the indeces of allData for that mouse's CNO data
            Sal{m} = currSal; %each entry corresponds to uniqueMice and contains the indeces of allData for that mouse's Saline data

            %Sort CNO condition into hours (regardless of day)
            [cnoDays, iFirst, iDays] = unique({allData(currCNO).date}); %find each day of CNO injection for this mouse
            hour1 = []; hour2 = []; hour3 = []; %will contain the file numbers for each hour for all days
            for d = 1:length(cnoDays) %for each day
                currDay = find(iDays == d); %indices of allData(currCNO) that correspond to current day
                currDay = currCNO(currDay); %indeces of alldata that correspond to current day's files
                currNums = {allData(currDay).tseriesNum};
                currNums = cellfun(@str2num, currNums, 'UniformOutput',0); %turn string into numbers
                currNums = cell2mat(currNums); %turn cell array into vector
                [~, i] = sort(currNums); %sort in order  

                currDaySorted = currDay(i); %sort file numbers in order of tSeries
                hours = NaN(1,3); %assume 3 hours for each day
                hours(1:length(currDay)) = currDaySorted; %if an hour wasn't recorded that day, will be NaN

                hour1 = [hour1 hours(1)]; %add file number to a vector (could be multiple CNO days, and so multiple CNO-hour1s per mouse)
                hour2 = [hour2 hours(2)];
                hour3 = [hour3 hours(3)];
            end
            CNO_1{m} = hour1; 
            CNO_2{m} = hour2; 
            CNO_3{m} = hour3;

            %Sort saline condition into hours 
            [salDays, iFirst, iDays] = unique({allData(currSal).date}); %find each day of CNO injection for this mouse
            hour1 = []; hour2 = []; hour3 = []; %will contain the file numbers for each hour for all days
            for d = 1:length(salDays) %for each day
                currDay = find(iDays == d); %indices of allData(currCNO) that correspond to current day
                currDay = currSal(currDay); %indeces of alldata that correspond to current day's files
                currNums = {allData(currDay).tseriesNum};
                currNums = cellfun(@str2num, currNums, 'UniformOutput',0); %turn string into numbers
                currNums = cell2mat(currNums); %turn cell array into vector
                [~, i] = sort(currNums); %sort in order  

                currDaySorted = currDay(i); %sort file numbers in order of tSeries
                hours = NaN(1,3); %assume 3 hours for each day
                hours(1:length(currDay)) = currDaySorted; %if an hour wasn't recorded that day, will be NaN

                hour1 = [hour1 hours(1)]; %add file number to a vector (could be multiple CNO days, and so multiple CNO-hour1s per mouse)
                hour2 = [hour2 hours(2)];
                hour3 = [hour3 hours(3)];
            end
            Sal_1{m} = hour1; 
            Sal_2{m} = hour2; 
            Sal_3{m} = hour3;
        end


        GiFiles = find(~cellfun(@isempty, strfind({allData(:).dreadd}, 'Gi')));
        GqFiles = find(~cellfun(@isempty, strfind({allData(:).dreadd}, 'Gq')));


        lfpMice = []; %corresponds to uniqueMice, which mice have lfp rec
        eegMice = []; %corresponds to uniqueMice, which mice have eeg rec
        for m = 1:length(uniqueMice)
            if ~isempty(allData(Sal_1{m}(1)).lfp) && ~isempty(allData(CNO_1{m}(1)).lfp)
                lfpMice = [lfpMice m];
            end
            if ~isempty(allData(Sal_1{m}(1)).eeg) && ~isempty(allData(CNO_1{m}(1)).eeg)
                eegMice = [eegMice m];
            end           
        end

        lfpFiles = [];
        eegFiles = [];
        for d = 1:length(allData)
            if ~isempty(allData(d).lfp)
                lfpFiles = [lfpFiles d];
            end
            if ~isempty(allData(d).eeg)
                eegFiles = [eegFiles d];
            end
        end    
    %% SLEEP-SCORE (Mult-files)
        sf = 1000; %samples per sec
        tic
        for m = 1:length(uniqueMice)
            files = mice{m};
            %sleep score based on LFP
            if ismember(m,lfpMice) %removed some lfp channels for noise, in that case use eeg
                data = {}; names = {};
                k=1;
                for f = Sal{m}'
                    temp = [allData(f).times, allData(f).lfp, allData(f).emg, round(allData(f).opto)];
                    data{k} = temp;
                    names{k} = allData(f).name;
                    k = k+1;
                end

                [INT, IDX, f_FFT, FFT, t_FFT, v_FFT, swTimes] = ...
                    sleepDetection_zscore_multFiles(data, sf, names, 0);

                k = 1;
                for f = Sal{m}'
                    allData(f).INT_lfp = (INT{k});
                    allData(f).IDX_lfp = single(IDX{k});
                    allData(f).f_FFT_lfp = single(f_FFT{k});
                    allData(f).FFT_lfp = single(FFT{k});
                    allData(f).t_FFT_lfp = single(t_FFT{k});
                    allData(f).v_FFT_lfp = single(v_FFT{k});
                    allData(f).swTimes_lfp = single(swTimes{k});
                    k=k+1;
                end

                data = {}; names = {};
                k = 1;
                for f = CNO{m}'
                    temp = [allData(f).times, allData(f).lfp, allData(f).emg, round(allData(f).opto)];
                    data{k} = temp;
                    names{k} = allData(f).name;
                    k = k+1;
                end

                [INT, IDX, f_FFT, FFT, t_FFT, v_FFT, swTimes] = ...
                    sleepDetection_zscore_multFiles(data, sf, names, 0);

                k = 1;
                for f = CNO{m}'
                    allData(f).INT_lfp = (INT{k});
                    allData(f).IDX_lfp = single(IDX{k});
                    allData(f).f_FFT_lfp = single(f_FFT{k});
                    allData(f).FFT_lfp = single(FFT{k});
                    allData(f).t_FFT_lfp = single(t_FFT{k});
                    allData(f).v_FFT_lfp = single(v_FFT{k});
                    allData(f).swTimes_lfp = single(swTimes{k});
                    k=k+1;
                end
            end  

            %sleep score based on EEG
            if ismember(m,eegMice) %if there is also an eeg channel
                data = {}; names = {};
                k=1;
                for f = Sal{m}'
                    temp = [allData(f).times, allData(f).eeg, allData(f).emg, round(allData(f).opto)];
                    data{k} = temp;
                    names{k} = allData(f).name;
                    k = k+1;
                end

                [INT, IDX, f_FFT, FFT, t_FFT, v_FFT, swTimes] = ...
                    sleepDetection_zscore_multFiles(data, sf, names, 0);

                k = 1;
                for f = Sal{m}'
                    allData(f).INT_eeg = (INT{k});
                    allData(f).IDX_eeg = single(IDX{k});
                    allData(f).f_FFT_eeg = single(f_FFT{k});
                    allData(f).FFT_eeg = single(FFT{k});
                    allData(f).t_FFT_eeg = single(t_FFT{k});
                    allData(f).v_FFT_eeg = single(v_FFT{k});
                    allData(f).swTimes_eeg = single(swTimes{k});
                    k=k+1;
                end

                data = {}; names = {};
                k=1;
                for f = CNO{m}'
                    temp = [allData(f).times, allData(f).eeg, allData(f).emg, round(allData(f).opto)];
                    data{k} = temp;
                    names{k} = allData(f).name;
                    k = k+1;
                end

                [INT, IDX, f_FFT, FFT, t_FFT, v_FFT, swTimes] = ...
                    sleepDetection_zscore_multFiles(data, sf, names, 0);

                k = 1;
                for f = CNO{m}'
                    allData(f).INT_eeg = (INT{k});
                    allData(f).IDX_eeg = single(IDX{k});
                    allData(f).f_FFT_eeg = single(f_FFT{k});
                    allData(f).FFT_eeg = single(FFT{k});
                    allData(f).t_FFT_eeg = single(t_FFT{k});
                    allData(f).v_FFT_eeg = single(v_FFT{k});
                    allData(f).swTimes_eeg = single(swTimes{k});
                    k=k+1;
                end
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

        %Iterate through freq bands and find average power over time using LFP
        for d = lfpFiles
            specSelect = []; %each row is power over time, corresponds to freqList
            for i = 1:size(freqList,1)
                first = freqList(i,1);
                last = freqList(i,2);

                [~,f_first] = min(abs(allData(d).f_FFT_lfp-first));
                [~,f_last] = min(abs(allData(d).f_FFT_lfp-last));

                temp = log(allData(d).FFT_lfp(:,f_first:f_last));
                specSelect(i,:) = mean(temp,2);

                %smooth
                specSelect(i,:) = smooth((specSelect(i,:)),20);
            end
            allData(d).specSelect_lfp = specSelect;
        end

        %Iterate through freq bands and find average power over time using EEG
        for d = eegFiles 
            specSelect = []; %each row is power over time, corresponds to freqList
            for i = 1:size(freqList,1)
                first = freqList(i,1);
                last = freqList(i,2);

                [~,f_first] = min(abs(allData(d).f_FFT_eeg-first));
                [~,f_last] = min(abs(allData(d).f_FFT_eeg-last));

                temp = log(allData(d).FFT_eeg(:,f_first:f_last));
                specSelect(i,:) = mean(temp,2);

                %smooth
                specSelect(i,:) = smooth((specSelect(i,:)),20);
            end
            allData(d).specSelect_eeg = specSelect;
        end
    %% Sleep score for baseline files
    %use sleep score to get spectral data
        sf = 1000; %samples per sec
        tic
        for d = 1:length(allData)
            if ~isempty(allData(d).dataBaseline)
                %sleep score based on LFP
                if ~isempty(allData(d).baseline_lfp) %removed some lfp channels for noise, in that case use eeg
                    temp = [allData(d).baseline_times, zscore(allData(d).baseline_lfp), allData(d).baseline_emg, allData(d).baseline_opto];

                    [~, ~, f_FFT, FFT, t_FFT, ~, ~] ...
                    = sleepDetection_zscore_051619(temp, sf, allData(d).name,0);

                    allData(d).baseline_f_FFT_lfp = single(f_FFT);
                    allData(d).baseline_FFT_lfp = single(FFT);
                    allData(d).baseline_t_FFT_lfp = single(t_FFT);
                end

                %sleep score based on EEG
                if ~isempty(allData(d).eeg) %if there is also an eeg channel
                    temp = [allData(d).baseline_times, zscore(allData(d).baseline_eeg), allData(d).baseline_emg, allData(d).baseline_opto];

                    [~, ~, f_FFT, FFT, t_FFT, ~, ~] ...
                    = sleepDetection_zscore_051619(temp, sf, allData(d).name,0);

                    allData(d).baseline_f_FFT_eeg = single(f_FFT);
                    allData(d).baseline_FFT_eeg = single(FFT);
                    allData(d).baseline_t_FFT_eeg = single(t_FFT);
                end
            end
            disp(strcat('done with file number',{' '},num2str(d)));
        end
        toc

        %find spectral info
        freqList = [.5 4 ;
        6 12 ;
        15 30 ;
        40 60 ;
        65 100];
        freqNames = {'delta', 'theta', 'beta', 'low gamma', 'high gamma'};

        %Iterate through freq bands and find average power over time using LFP
        for d = 1:length(allData)
            if ~isempty(allData(d).dataBaseline)
                specSelectLFP = []; %each row is power over time, corresponds to freqList
                specSelectEEG = [];
                for i = 1:size(freqList,1)            
                    if ~isempty(allData(d).baseline_lfp)
                        first = freqList(i,1);
                        last = freqList(i,2);

                        [~,f_first] = min(abs(allData(d).baseline_f_FFT_lfp-first));
                        [~,f_last] = min(abs(allData(d).baseline_f_FFT_lfp-last));

                        temp = log(allData(d).baseline_FFT_lfp(:,f_first:f_last));
                        specSelectLFP(i,:) = mean(temp,2);

                        %smooth
                        specSelectLFP(i,:) = smooth((specSelectLFP(i,:)),20);
                    end

                    if ~isempty(allData(d).baseline_eeg)
                        first = freqList(i,1);
                        last = freqList(i,2);

                        [~,f_first] = min(abs(allData(d).baseline_f_FFT_eeg-first));
                        [~,f_last] = min(abs(allData(d).baseline_f_FFT_eeg-last));

                        temp = log(allData(d).baseline_FFT_eeg(:,f_first:f_last));
                        specSelectEEG(i,:) = mean(temp,2);

                        %smooth
                        specSelectEEG(i,:) = smooth((specSelectEEG(i,:)),20);
                    end
                allData(d).baseline_specSelect_lfp = specSelectLFP;
                allData(d).baseline_specSelect_eeg = specSelectEEG;
                end
            end
        end

        hasbaseline = [];
        for m = 1:length(uniqueMice)
            currSal = Sal_1{m}; currCNO = CNO_1{m};
            if ~isempty(allData(currSal).dataBaseline) && ~isempty(allData(currCNO).dataBaseline)
                hasbaseline = [hasbaseline m];
            end
        end   
    %% Calculate locomotion
    %%%%% Movement buffer 5 s here, but 10 s for endgoenous dataset.
        win = 1000; %window in ms to look at opto
        movementBufferWin = 5; %in seconds: how many sec before start of stationary period must there be no locomotion

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

            if ~isempty(allData(d).INT_lfp); wakes = allData(d).INT_lfp{1,1}; else; wakes = allData(d).INT_eeg{1,1}; end
            for i = 1:size(wakes,1)
                start = wakes(i,1); %when wake period starts in sec
                finish = wakes(i,2); %finish of wake period in sec
                start = round(start);
                finish = round(finish);

                temp = find(allData(d).locomotion(start:finish)) + (start-1); %finds indices of locomotion where locomotion is 1 (ie the sec where there is movement)
                if ~isempty(temp)
                    move = find(diff(temp)>1); 
                    e = [temp(move), temp(end)]; %marks all the end points of the movement period
                    s = [temp(1) temp(move + 1)]; %marks all the start points of the movement period

                    allData(d).whenMoving = [allData(d).whenMoving ; s' e'];
                end

                temp = find(allData(d).locomotion(start:finish) == 0) + (start-1); %finds idx where locomotion is 0 (ie the sec where there is movement)
                if ~isempty(temp)
                    Nomove = find(diff(temp)>1); %when there are >1 period 
                    e = [temp(Nomove), temp(end)]; %marks all the end points of the movement period
                    s = [temp(1) temp(Nomove + 1)]; %marks all the start points of the movement period

            %Add condition here where non-moving periods cannot immediately preceed a moving period - must have a buffer (movementBufferWin)        
                    lengths = e-s;
                    toKeep = lengths > 10;

                    s = s(toKeep);
                    e = e(toKeep);

                    s = s + movementBufferWin; 

                    allData(d).whenNotMoving = [allData(d).whenNotMoving ; s' e'];
                end
            end
        end
    %% Calculate event rate
        for d = 1:length(allRes)
            if ~isempty(allRes(d).opts)
                allRes(d).frameperiod = allRes(d).opts.frameRate;
                allRes(d).numFrames = allRes(d).opts.sz(3);
                allRes(d).peakTimesF = allRes(d).fts.loc.t0; %in frames (onsets)
                allRes(d).peakTimes = allRes(d).peakTimesF*allRes(d).frameperiod; %in sec
                allData(d).peakTimes = allRes(d).peakTimes;
                [allRes(d).numEvents, allRes(d).timebins] = hist(allRes(d).peakTimes, 1:(allRes(d).numFrames * allRes(d).frameperiod)); %
                allRes(d).eventRate = allRes(d).numEvents/(allRes(d).timebins(2)-allRes(d).timebins(1)); %events per second
            end
        end
        for d = 1:length(allData)
            [allData(d).onsetSec, allData(d).EvtsSortByOnset] = sort(allData(d).peakTimes); %EvtsSortByOnset is list of event IDs in order of when they happen
            allRes(d).onsetSec = allData(d).onsetSec;
            allRes(d).EvtsSortByOnset = allData(d).EvtsSortByOnset;
        end
    %% Analyze based on behavioral state: sleep    

    minSWS = 15; %in sec: min duration of sleep period to be used
    for d = 1:length(allData)
        if ~isempty(allRes(d).opts)
            if ~isempty(allData(d).INT_lfp)
                sleepNumEvents = [];
                sleepEventRate = [];
                sleepDffMax = [];
                sleepDffMax_PeriodMean = []; %each entry is the mean dff for that sleep period
                sleepArea = []; %area in pixels
                sleepArea_PeriodMean = []; %each entry is the mean size (in pixels) for that sleep period
                sleepDuration = []; %duration of events in frames (peak 10% to peak 10%)
                sleepDuration_PeriodMean = []; %each entry is mean duration (in frames) for that sleep period

                sleepDelta_PeriodMean = []; %each entry is mean delta power for that sleep period
                sleepTheta_PeriodMean = []; %each entry is mean theta power for that sleep period

                allData(d).sleepEvtIdx_lfp = {};

                for i = 1:length(allData(d).INT_lfp{1,2})
                    start = allData(d).INT_lfp{1,2}(i,1); %when sleep period starts in sec
                    finish = allData(d).INT_lfp{1,2}(i,2); %when sleep period ends in sec

                    if finish-start > minSWS && finish<length(allRes(d).numEvents)

                        startIdx = start; %find idx of timebins corresponding to start of sleep period
                        finishIdx = finish; %find idx of timebins corresponding to end of sleep period

                        sleepNumEvents = [sleepNumEvents sum(allRes(d).numEvents(startIdx:finishIdx))];
                        sleepEventRate = [sleepEventRate mean(allRes(d).eventRate(startIdx:finishIdx))];

                        relEvents = allRes(d).EvtsSortByOnset(find(allRes(d).onsetSec >= start & allRes(d).onsetSec <= finish)); %event IDs that occur in this period
                        allData(d).sleepEvtIdx_lfp{i} = relEvents; 

                        sleepDffMax = [sleepDffMax allRes(d).fts.curve.dffMax(relEvents)];
                        sleepDffMax_PeriodMean = [sleepDffMax_PeriodMean mean(allRes(d).fts.curve.dffMax(relEvents))];
                        sleepArea = [sleepArea allRes(d).fts.basic.area(relEvents)];
                        sleepArea_PeriodMean = [sleepArea_PeriodMean mean(allRes(d).fts.basic.area(relEvents))];
                        sleepDuration = [sleepDuration allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents)];
                        sleepDuration_PeriodMean = [sleepDuration_PeriodMean mean(allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents))];
                    end
                end

                allData(d).sleepNumEvents_lfp = sleepNumEvents;
                allData(d).sleepEventRate_lfp = sleepEventRate;
                allData(d).sleepDffMax_lfp = sleepDffMax;
                allData(d).sleepDffMax_PeriodMean_lfp = sleepDffMax_PeriodMean; %each entry is the mean dff for that sleep period
                allData(d).sleepArea_lfp = sleepArea; %area in pixels
                allData(d).sleepArea_PeriodMean_lfp = sleepArea_PeriodMean; %each entry is the mean size (in pixels) for that sleep period
                allData(d).sleepDuration_lfp = sleepDuration; %duration of events in frames (peak 10% to peak 10%)
                allData(d).sleepDuration_PeriodMean_lfp = sleepDuration_PeriodMean; %each entry is mean duration (in frames) for that sleep period
            end

            if ~isempty(allData(d).INT_eeg)
                sleepNumEvents = [];
                sleepEventRate = [];
                sleepDffMax = [];
                sleepDffMax_PeriodMean = []; %each entry is the mean dff for that sleep period
                sleepArea = []; %area in pixels
                sleepArea_PeriodMean = []; %each entry is the mean size (in pixels) for that sleep period
                sleepDuration = []; %duration of events in frames (peak 10% to peak 10%)
                sleepDuration_PeriodMean = []; %each entry is mean duration (in frames) for that sleep period

                sleepDelta_PeriodMean = []; %each entry is mean delta power for that sleep period
                sleepTheta_PeriodMean = []; %each entry is mean theta power for that sleep period

                allData(d).sleepEvtIdx_eeg = {};

                for i = 1:length(allData(d).INT_eeg{1,2})
                    start = allData(d).INT_eeg{1,2}(i,1); %when sleep period starts in sec
                    finish = allData(d).INT_eeg{1,2}(i,2); %when sleep period ends in sec

                    if finish-start > minSWS && finish<length(allRes(d).numEvents)

                        startIdx = start; %find idx of timebins corresponding to start of sleep period
                        finishIdx = finish; %find idx of timebins corresponding to end of sleep period

                        sleepNumEvents = [sleepNumEvents sum(allRes(d).numEvents(startIdx:finishIdx))];
                        sleepEventRate = [sleepEventRate mean(allRes(d).eventRate(startIdx:finishIdx))];

                        relEvents = allRes(d).EvtsSortByOnset(find(allRes(d).onsetSec >= start & allRes(d).onsetSec <= finish)); %event IDs that occur in this period
                        allData(d).sleepEvtIdx_eeg{i} = relEvents; 

                        sleepDffMax = [sleepDffMax allRes(d).fts.curve.dffMax(relEvents)];
                        sleepDffMax_PeriodMean = [sleepDffMax_PeriodMean mean(allRes(d).fts.curve.dffMax(relEvents))];
                        sleepArea = [sleepArea allRes(d).fts.basic.area(relEvents)];
                        sleepArea_PeriodMean = [sleepArea_PeriodMean mean(allRes(d).fts.basic.area(relEvents))];
                        sleepDuration = [sleepDuration allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents)];
                        sleepDuration_PeriodMean = [sleepDuration_PeriodMean mean(allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents))];
                    end
                end

                allData(d).sleepNumEvents_eeg = sleepNumEvents;
                allData(d).sleepEventRate_eeg = sleepEventRate;
                allData(d).sleepDffMax_eeg = sleepDffMax;
                allData(d).sleepDffMax_PeriodMean_eeg = sleepDffMax_PeriodMean; %each entry is the mean dff for that sleep period
                allData(d).sleepArea_eeg = sleepArea; %area in pixels
                allData(d).sleepArea_PeriodMean_eeg = sleepArea_PeriodMean; %each entry is the mean size (in pixels) for that sleep period
                allData(d).sleepDuration_eeg = sleepDuration; %duration of events in frames (peak 10% to peak 10%)
                allData(d).sleepDuration_PeriodMean_eeg = sleepDuration_PeriodMean; %each entry is mean duration (in frames) for that sleep period
            end
        end
    end
    %% Analyze based on behavioral state: wake   

    for d = 1:length(allData)
        if ~isempty(allRes(d).opts)
            if ~isempty(allData(d).INT_lfp)
                wakeNumEvents = [];
                wakeEventRate = [];
                wakeDffMax = [];
                wakeDffMax_PeriodMean = []; %each entry is the mean dff for that wake period
                wakeArea = []; %area in pixels
                wakeArea_PeriodMean = []; %each entry is the mean size (in pixels) for that wake period
                wakeDuration = []; %duration of events in frames (peak 10% to peak 10%)
                wakeDuration_PeriodMean = []; %each entry is mean duration (in frames) for that wake period

                wakeDelta_PeriodMean = []; %each entry is mean delta power for that sleep period
                wakeTheta_PeriodMean = []; %each entry is mean theta power for that sleep period

                allData(d).wakeEvtIdx_lfp = {};

                for i = 1:length(allData(d).INT_lfp{1,1})
                    start = allData(d).INT_lfp{1,1}(i,1); %when sleep period starts in sec
                    finish = allData(d).INT_lfp{1,1}(i,2); %when sleep period ends in sec

                    startIdx = start; %find idx of timebins corresponding to start of sleep period
                    finishIdx = finish; %find idx of timebins corresponding to end of sleep period

                    if finish<length(allRes(d).numEvents)
                        wakeNumEvents = [wakeNumEvents sum(allRes(d).numEvents(startIdx:finishIdx))];
                        wakeEventRate = [wakeEventRate mean(allRes(d).eventRate(startIdx:finishIdx))];

                        relEvents = allRes(d).EvtsSortByOnset(allRes(d).onsetSec >= start & allRes(d).onsetSec <= finish); %event IDs that occur in this period
                        allData(d).wakeEvtIdx_lfp{i} = relEvents;

                        wakeDffMax = [wakeDffMax allRes(d).fts.curve.dffMax(relEvents)];
                        wakeDffMax_PeriodMean = [wakeDffMax_PeriodMean mean(allRes(d).fts.curve.dffMax(relEvents))];
                        wakeArea = [wakeArea allRes(d).fts.basic.area(relEvents)];
                        wakeArea_PeriodMean = [wakeArea_PeriodMean mean(allRes(d).fts.basic.area(relEvents))];
                        wakeDuration = [wakeDuration allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents)];
                        wakeDuration_PeriodMean = [wakeDuration_PeriodMean mean(allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents))];
                    end
                end

                allData(d).wakeNumEvents_lfp = wakeNumEvents;
                allData(d).wakeEventRate_lfp = wakeEventRate;
                allData(d).wakeDffMax_lfp = wakeDffMax;
                allData(d).wakeDffMax_PeriodMean_lfp = wakeDffMax_PeriodMean; %each entry is the mean dff for that wake period
                allData(d).wakeArea_lfp = wakeArea; %area in pixels
                allData(d).wakeArea_PeriodMean_lfp = wakeArea_PeriodMean; %each entry is the mean size (in pixels) for that wake period
                allData(d).wakeDuration_lfp = wakeDuration; %duration of events in frames (peak 10% to peak 10%)
                allData(d).wakeDuration_PeriodMean_lfp = wakeDuration_PeriodMean; %each entry is mean duration (in frames) for that wake period
            end

            if ~isempty(allData(d).INT_eeg)
                wakeNumEvents = [];
                wakeEventRate = [];
                wakeDffMax = [];
                wakeDffMax_PeriodMean = []; %each entry is the mean dff for that wake period
                wakeArea = []; %area in pixels
                wakeArea_PeriodMean = []; %each entry is the mean size (in pixels) for that wake period
                wakeDuration = []; %duration of events in frames (peak 10% to peak 10%)
                wakeDuration_PeriodMean = []; %each entry is mean duration (in frames) for that wake period

                wakeDelta_PeriodMean = []; %each entry is mean delta power for that sleep period
                wakeTheta_PeriodMean = []; %each entry is mean theta power for that sleep period

                allData(d).wakeEvtIdx_eeg = {};

                for i = 1:length(allData(d).INT_eeg{1,1})
                    start = allData(d).INT_eeg{1,1}(i,1); %when sleep period starts in sec
                    finish = allData(d).INT_eeg{1,1}(i,2); %when sleep period ends in sec

                    startIdx = start; %find idx of timebins corresponding to start of sleep period
                    finishIdx = finish; %find idx of timebins corresponding to end of sleep period

                    if finish<length(allRes(d).numEvents)

                        wakeNumEvents = [wakeNumEvents sum(allRes(d).numEvents(startIdx:finishIdx))];
                        wakeEventRate = [wakeEventRate mean(allRes(d).eventRate(startIdx:finishIdx))];

                        relEvents = allRes(d).EvtsSortByOnset(allRes(d).onsetSec >= start & allRes(d).onsetSec <= finish); %event IDs that occur in this period
                        allData(d).wakeEvtIdx_eeg{i} = relEvents;

                        wakeDffMax = [wakeDffMax allRes(d).fts.curve.dffMax(relEvents)];
                        wakeDffMax_PeriodMean = [wakeDffMax_PeriodMean mean(allRes(d).fts.curve.dffMax(relEvents))];
                        wakeArea = [wakeArea allRes(d).fts.basic.area(relEvents)];
                        wakeArea_PeriodMean = [wakeArea_PeriodMean mean(allRes(d).fts.basic.area(relEvents))];
                        wakeDuration = [wakeDuration allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents)];
                        wakeDuration_PeriodMean = [wakeDuration_PeriodMean mean(allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents))];
                    end
                end

                allData(d).wakeNumEvents_eeg = wakeNumEvents;
                allData(d).wakeEventRate_eeg = wakeEventRate;
                allData(d).wakeDffMax_eeg = wakeDffMax;
                allData(d).wakeDffMax_PeriodMean_eeg = wakeDffMax_PeriodMean; %each entry is the mean dff for that wake period
                allData(d).wakeArea_eeg = wakeArea; %area in pixels
                allData(d).wakeArea_PeriodMean_eeg = wakeArea_PeriodMean; %each entry is the mean size (in pixels) for that wake period
                allData(d).wakeDuration_eeg = wakeDuration; %duration of events in frames (peak 10% to peak 10%)
                allData(d).wakeDuration_PeriodMean_eeg = wakeDuration_PeriodMean; %each entry is mean duration (in frames) for that wake period
            end        

        end
    end
    %% Analyze based on behavioral state: wake-moving 

    for d = 1:length(allData)
        if ~isempty(allRes(d).opts)
            movingNumEvents = [];
            movingEventRate = [];
            movingDffMax = [];
            movingDffMax_PeriodMean = []; %each entry is the mean dff for that wake period
            movingArea = []; %area in pixels
            movingArea_PeriodMean = []; %each entry is the mean size (in pixels) for that wake period
            movingDuration = []; %duration of events in frames (peak 10% to peak 10%)
            movingDuration_PeriodMean = []; %each entry is mean duration (in frames) for that wake period

            movingDelta_PeriodMean = []; %each entry is mean delta power for that sleep period
            movingTheta_PeriodMean = []; %each entry is mean theta power for that sleep period

            allData(d).movingEvtIdx = {}; 

            for i = 1:length(allData(d).whenMoving)
                start = allData(d).whenMoving(i,1); %when moving period starts in sec
                finish = allData(d).whenMoving(i,2); %when moving period ends in sec

                startIdx = start; %find idx of timebins corresponding to start of sleep period
                finishIdx = finish; %find idx of timebins corresponding to end of sleep period

                if finish<length(allRes(d).numEvents)
                    movingNumEvents = [movingNumEvents sum(allRes(d).numEvents(startIdx:finishIdx))];
                    movingEventRate = [movingEventRate mean(allRes(d).eventRate(startIdx:finishIdx))];

                    relEvents = allRes(d).EvtsSortByOnset(find(allRes(d).onsetSec >= start & allRes(d).onsetSec <= finish)); %event IDs that occur in this period
                    allData(d).movingEvtIdx{i} = relEvents;

                    movingDffMax = [movingDffMax allRes(d).fts.curve.dffMax(relEvents)];
                    movingDffMax_PeriodMean = [movingDffMax_PeriodMean mean(allRes(d).fts.curve.dffMax(relEvents))];
                    movingArea = [movingArea allRes(d).fts.basic.area(relEvents)];
                    movingArea_PeriodMean = [movingArea_PeriodMean mean(allRes(d).fts.basic.area(relEvents))];
                    movingDuration = [movingDuration allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents)];
                    movingDuration_PeriodMean = [movingDuration_PeriodMean mean(allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents))];
                end
            end

            allData(d).movingNumEvents = movingNumEvents;
            allData(d).movingEventRate = movingEventRate;
            allData(d).movingDffMax = movingDffMax;
            allData(d).movingDffMax_PeriodMean = movingDffMax_PeriodMean; %each entry is the mean dff for that wake period
            allData(d).movingArea = movingArea; %area in pixels
            allData(d).movingArea_PeriodMean = movingArea_PeriodMean; %each entry is the mean size (in pixels) for that wake period
            allData(d).movingDuration = movingDuration; %duration of events in frames (peak 10% to peak 10%)
            allData(d).movingDuration_PeriodMean = movingDuration_PeriodMean; %each entry is mean duration (in frames) for that wake period
        end
    end
    %% Analyze based on behavioral state: wake-stationary

    minStationary = 15; %in sec: min duration of a stationary period
    for d = 1:length(allData)
        if ~isempty(allRes(d).opts)
            notmovingNumEvents = [];
            notmovingEventRate = [];
            notmovingDffMax = [];
            notmovingDffMax_PeriodMean = []; %each entry is the mean dff for that wake period
            notmovingArea = []; %area in pixels
            notmovingArea_PeriodMean = []; %each entry is the mean size (in pixels) for that wake period
            notmovingDuration = []; %duration of events in frames (peak 10% to peak 10%)
            notmovingDuration_PeriodMean = []; %each entry is mean duration (in frames) for that wake period

            notmovingDelta_PeriodMean = []; %each entry is mean delta power for that sleep period
            notmovingTheta_PeriodMean = []; %each entry is mean theta power for that sleep period

            allData(d).notmovingEvtIdx = {};

            for i = 1:length(allData(d).whenNotMoving)
                start = allData(d).whenNotMoving(i,1); %when stationary period starts in sec
                finish = allData(d).whenNotMoving(i,2); %when stationary period ends in sec

                if finish-start > minStationary

                    startIdx = start; %find idx of timebins corresponding to start of sleep period
                    finishIdx = finish; %find idx of timebins corresponding to end of sleep period

                    if finish<length(allRes(d).numEvents)
                        notmovingNumEvents = [notmovingNumEvents sum(allRes(d).numEvents(startIdx:finishIdx))];
                        notmovingEventRate = [notmovingEventRate mean(allRes(d).eventRate(startIdx:finishIdx))];

                        relEvents = allRes(d).EvtsSortByOnset(find(allRes(d).onsetSec >= start & allRes(d).onsetSec <= finish)); %event IDs that occur in this period

                        allData(d).notmovingEvtIdx{i} = relEvents;

                        notmovingDffMax = [notmovingDffMax allRes(d).fts.curve.dffMax(relEvents)];
                        notmovingDffMax_PeriodMean = [notmovingDffMax_PeriodMean mean(allRes(d).fts.curve.dffMax(relEvents))];
                        notmovingArea = [notmovingArea allRes(d).fts.basic.area(relEvents)];
                        notmovingArea_PeriodMean = [notmovingArea_PeriodMean mean(allRes(d).fts.basic.area(relEvents))];
                        notmovingDuration = [notmovingDuration allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents)];
                        notmovingDuration_PeriodMean = [notmovingDuration_PeriodMean mean(allRes(d).fts.curve.tEnd(relEvents)-allRes(d).fts.curve.tBegin(relEvents))];
                    end
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
        end
    end
    %% Slow oscillations vs delta
    %Find SOs and deltas for all files using detectSOvDelta_121319
    
    %First create a vector, sleepTimes, with 0/1 corresponding to sleep for
    %every second
    for d = 1:length(allData)
        recLength = round(allData(d).times(end)/1000);
        if ~isempty(allData(d).INT_lfp)
            sleeps = round(allData(d).INT_lfp{2});
            sleepTimes = zeros(1, recLength); 
            for s = 1:size(sleeps,1)
                sleepTimes(sleeps(s,1):sleeps(s,2)) = 1;
            end
            allData(d).sleepTimes_lfp = sleepTimes;
        end
    end
    
    
    for d = 1:length(allData)
        if ~isempty(allData(d).sleepTimes_lfp)
            [out, lfp_delta, lfp_gamma, flip] = detectSOvDelta_121319(double(zscore((allData(d).lfp))),1000,allData(d).sleepTimes_lfp,0);
            allData(d).SOstruct_lfp = out;
            allData(d).lfp_gamma = lfp_gamma;
            allData(d).lfp_delta = lfp_delta;
            if flip
                allData(d).lfp = allData(d).lfp * -1;
                allData(d).lfp_gamma = allData(d).lfp_gamma * -1; 
            end
            disp(strcat('done with LFP file',{' '},num2str(d),{' '},'flip =',num2str(flip)))
        end
    end    

%% Fig4FigSupp1 A LEFT, C LEFT (peak CDFs)

minDur = 10;

allPeaksD_sal = []; allPeaks_sal = []; %each entry is the mean value of a sleep period
allTroughsD_sal = []; allTroughs_sal = [];
for d = [Sal_1{:} Sal_2{:}]
    if ~isempty(allData(d).SOstruct_lfp)
        sleeps = allData(d).INT_lfp{2};
        sLengths = sleeps(:,2)-sleeps(:,1);
        sleeps = sleeps(sLengths>minDur,:);
        for s = 1:size(sleeps,1)
            currD = find(allData(d).SOstruct_lfp.zc_d >= sleeps(s,1) & ...
                allData(d).SOstruct_lfp.zc_d < sleeps(s,2));
            currSO = find(allData(d).SOstruct_lfp.zc >= sleeps(s,1) & ...
                allData(d).SOstruct_lfp.zc < sleeps(s,2));

            allPeaksD_sal = [allPeaksD_sal ; mean(allData(d).SOstruct_lfp.peaks_d(currD))];
            allPeaks_sal = [allPeaks_sal ; mean(allData(d).SOstruct_lfp.peaks(currSO))];
            allTroughsD_sal = [allTroughsD_sal ; mean(allData(d).SOstruct_lfp.troughs_d(currD))];
            allTroughs_sal = [allTroughs_sal ; mean(allData(d).SOstruct_lfp.troughs(currSO))];
        end
    end
end

allPeaksD_cno = []; allPeaks_cno = [];
allTroughsD_cno = []; allTroughs_cno = [];
for d = [CNO_1{:} CNO_2{:}]
    if ~isempty(allData(d).SOstruct_lfp)
        sleeps = allData(d).INT_lfp{2};
        sLengths = sleeps(:,2)-sleeps(:,1);
        sleeps = sleeps(sLengths>minDur,:);
        for s = 1:size(sleeps,1)
            currD = find(allData(d).SOstruct_lfp.zc_d >= sleeps(s,1) & ...
                allData(d).SOstruct_lfp.zc_d < sleeps(s,2));
            currSO = find(allData(d).SOstruct_lfp.zc >= sleeps(s,1) & ...
                allData(d).SOstruct_lfp.zc < sleeps(s,2));

            allPeaksD_cno = [allPeaksD_cno ; mean(allData(d).SOstruct_lfp.peaks_d(currD))];
            allPeaks_cno = [allPeaks_cno ; mean(allData(d).SOstruct_lfp.peaks(currSO))];
            allTroughsD_cno = [allTroughsD_cno ; mean(allData(d).SOstruct_lfp.troughs_d(currD))];
            allTroughs_cno = [allTroughs_cno ; mean(allData(d).SOstruct_lfp.troughs(currSO))];
        end
    end
end

%%% Figure Fig4FigSupp1A LEFT and Fig4FigSupp1C LEFT (peak):
[salVals, salXs] = histcounts(allPeaksD_sal,80,'Normalization','cdf');
[cnoVals, cnoXs] = histcounts(allPeaksD_cno,80,'Normalization','cdf');
figure;
subplot(1,2,1) %PANEL A LEFT
plot(salXs(1:end-1), salVals,'k')
hold on;
plot(cnoXs(1:end-1), cnoVals,'r')
title('delta'); ylabel('cumulative probability'); xlabel('peak amplitude')
[h,p,stat] = kstest2(allPeaksD_sal,allPeaksD_cno);
text(mean(salXs),.7,strcat('kstest2 p =',{' '},num2str(p)))
set(gca,'fontsize',15)

[salVals, salXs] = histcounts(allPeaks_sal,80,'Normalization','cdf');
[cnoVals, cnoXs] = histcounts(allPeaks_cno,80,'Normalization','cdf');
subplot(1,2,2); %PANEL C LEFT
plot(salXs(1:end-1), salVals,'k')
hold on;
plot(cnoXs(1:end-1), cnoVals,'r')
title('slow oscillations'); ylabel('cumulative probability'); xlabel('peak amplitude')
[h,p,stat] = kstest2(allPeaks_sal,allPeaks_cno);
text(mean(salXs),.7,strcat('kstest2 p =',{' '},num2str(p)))
set(gca,'fontsize',15)

%% Fig4FigSupp1 B LEFT, D LEFT (trough CDFs)

%%% Using same variables created above (allTroughs_ ...)

[salVals, salXs] = histcounts(allTroughsD_sal,80,'Normalization','cdf');
[cnoVals, cnoXs] = histcounts(allTroughsD_cno,80,'Normalization','cdf');
figure;
subplot(1,2,1) %PANEL B LEFT
plot(salXs(1:end-1), salVals,'k')
hold on;
plot(cnoXs(1:end-1), cnoVals,'r')
title('delta'); ylabel('cumulative probability'); xlabel('trough amplitude')
[h,p,stat] = kstest2(allTroughsD_sal,allTroughsD_cno);
text(mean(salXs),.7,strcat('kstest2 p =',{' '},num2str(p)))
set(gca,'fontsize',15)

[salVals, salXs] = histcounts(allTroughs_sal,80,'Normalization','cdf');
[cnoVals, cnoXs] = histcounts(allTroughs_cno,80,'Normalization','cdf');
subplot(1,2,2); %PANEL D LEFT
plot(salXs(1:end-1), salVals,'k')
hold on;
plot(cnoXs(1:end-1), cnoVals,'r')
title('slow oscillations'); ylabel('cumulative probability'); xlabel('trough amplitude')
[h,p,stat] = kstest2(allTroughs_sal,allTroughs_cno);
text(mean(salXs),.7,strcat('kstest2 p =',{' '},num2str(p)))
set(gca,'fontsize',15)

%% Fig4FigSupp1 A RIGHT, C RIGHT (peak errorbars)

miceToUse = lfpMice; 
deltaCNO = {}; deltaSal = {};
soCNO = {}; soSal = {};
for m = miceToUse
    downDeltaCNO = []; downDeltaSal = []; 
    downSOCNO = []; downSOSal = []; 
    for f = CNO{m}'
        downDeltaCNO = [downDeltaCNO ; mean(allData(f).SOstruct_lfp.peaks_d)];
        downSOCNO = [downSOCNO ; mean(allData(f).SOstruct_lfp.peaks)];
    end
    deltaCNO{m} = downDeltaCNO;
    soCNO{m} = downSOCNO;
    for f = Sal{m}'
        downDeltaSal = [downDeltaSal ; mean(allData(f).SOstruct_lfp.peaks_d)];
        downSOSal = [downSOSal ; mean(allData(f).SOstruct_lfp.peaks)];
    end
    deltaSal{m} = downDeltaSal;
    soSal{m} = downSOSal;
end

colors = jet(length(miceToUse));
figure;
subplot(1,2,1); %PANEL A RIGHT
    k=1;
    clear s
    for m = miceToUse
        s(k)=scatter([1 2], [mean(deltaSal{m}) mean(deltaCNO{m})],80,colors(k,:),'filled');
        hold on;
        errorbar([1 2], [mean(deltaSal{m}) mean(deltaCNO{m})], [0 0],'k');
        hold on;
        k=k+1;
    end
    sal = cellfun(@mean, deltaSal(miceToUse));
    cno = cellfun(@mean, deltaCNO(miceToUse));
    scatter([0.75 2.25], [mean(sal) mean(cno)],'k','filled')
    hold on;
    errorbar([0.75 2.25], [mean(sal) mean(cno)],[std(sal) std(cno)],...
        'k','linestyle','none')
    [h,p] = ttest(cellfun(@mean,deltaSal(miceToUse)), cellfun(@mean, deltaCNO(miceToUse)));
    hold on;
    text(1.2, max([sal cno]), strcat('p =',{' '},num2str(p)))
    xlim([0 3]); 
    xticks([1 2])
    xticklabels({'saline','CNO'}); 
    ylabel('mean amplitude of delta down state');
    title('Gi: LFP delta'); 
    set(gca,'fontsize',15)
subplot(1,2,2); %PANEL C RIGHT
    k=1;
    clear s
    for m = miceToUse
        s(k)=scatter([1 2], [mean(soSal{m}) mean(soCNO{m})],80,colors(k,:),'filled');
        hold on;
        errorbar([1 2], [mean(soSal{m}) mean(soCNO{m})], [0 0],'k');
        hold on;
        k=k+1;
    end
    sal = cellfun(@mean, soSal(miceToUse));
    cno = cellfun(@mean, soCNO(miceToUse));
    scatter([0.75 2.25], [mean(sal) mean(cno)],'k','filled')
    hold on;
    errorbar([0.75 2.25], [mean(sal) mean(cno)],[std(sal) std(cno)],...
        'k','linestyle','none')
    [h,p] = ttest(cellfun(@mean,soSal(miceToUse)), cellfun(@mean, soCNO(miceToUse)));
    hold on;
    text(1.2, max([sal cno]), strcat('p =',{' '},num2str(p)))
    xlim([0 3]); 
    xticks([1 2])
    xticklabels({'saline','CNO'}); 
    ylabel('mean amplitude of SO down state');
    title('Gi: LFP slow oscillation');
    set(gca,'fontsize',15)
    set(gcf,'renderer','painters')
    
%% Fig4FigSupp1 B RIGHT, D RIGHT (trough errorbars)

miceToUse = lfpMice; 
deltaCNO = {}; deltaSal = {};
soCNO = {}; soSal = {};
for m = miceToUse
    downDeltaCNO = []; downDeltaSal = []; %delta down state density normalized to sec of sleep
    downSOCNO = []; downSOSal = []; %SO down state density normalized to sec of sleep
    for f = CNO{m}'
        downDeltaCNO = [downDeltaCNO ; mean(allData(f).SOstruct_lfp.troughs_d)];
        downSOCNO = [downSOCNO ; mean(allData(f).SOstruct_lfp.troughs)];
    end
    deltaCNO{m} = downDeltaCNO;
    soCNO{m} = downSOCNO;
    for f = Sal{m}'
        downDeltaSal = [downDeltaSal ; mean(allData(f).SOstruct_lfp.troughs_d)];
        downSOSal = [downSOSal ; mean(allData(f).SOstruct_lfp.troughs)];
    end
    deltaSal{m} = downDeltaSal;
    soSal{m} = downSOSal;
end

colors = jet(length(miceToUse));
figure;
subplot(1,2,1); %PANEL B RIGHT
    k=1;
    clear s
    for m = miceToUse
        s(k)=scatter([1 2], [mean(deltaSal{m}) mean(deltaCNO{m})],80,colors(k,:),'filled');
        hold on;
        errorbar([1 2], [mean(deltaSal{m}) mean(deltaCNO{m})], [0 0],'k');
        hold on;
        k=k+1;
    end
    sal = cellfun(@mean, deltaSal(miceToUse));
    cno = cellfun(@mean, deltaCNO(miceToUse));
    scatter([0.75 2.25], [mean(sal) mean(cno)],'k','filled')
    hold on;
    errorbar([0.75 2.25], [mean(sal) mean(cno)],[std(sal) std(cno)],...
        'k','linestyle','none')
    [h,p] = ttest(cellfun(@mean,deltaSal(miceToUse)), cellfun(@mean, deltaCNO(miceToUse)));
    hold on;
    text(1.2, max([sal cno]), strcat('p =',{' '},num2str(p)))
    xlim([0 3]); 
    xticks([1 2])
    xticklabels({'saline','CNO'}); 
    ylabel('mean amplitude of delta up state');
    title('Gi: LFP delta');
    set(gca,'fontsize',15);ylim([-1.6 -1])
subplot(1,2,2); %PANEL D RIGHT
    k=1;
    clear s
    for m = miceToUse
        s(k)=scatter([1 2], [mean(soSal{m}) mean(soCNO{m})],80,colors(k,:),'filled');
        hold on;
        errorbar([1 2], [mean(soSal{m}) mean(soCNO{m})], [0 0],'k');
        hold on;
        k=k+1;
    end
    sal = cellfun(@mean, soSal(miceToUse));
    cno = cellfun(@mean, soCNO(miceToUse));
    scatter([0.75 2.25], [mean(sal) mean(cno)],'k','filled')
    hold on;
    errorbar([0.75 2.25], [mean(sal) mean(cno)],[std(sal) std(cno)],...
        'k','linestyle','none')
    [h,p] = ttest(cellfun(@mean,soSal(miceToUse)), cellfun(@mean, soCNO(miceToUse)));
    hold on;
    text(1.2, max([sal cno]), strcat('p =',{' '},num2str(p)))
    xlim([0 3]); 
    xticks([1 2])
    xticklabels({'saline','CNO'}); 
    ylabel('mean amplitude of SO up state');
    title('Gi: LFP slow oscillation');
    set(gca,'fontsize',15); ylim([-2 -1])
    set(gcf,'renderer','painters')
    