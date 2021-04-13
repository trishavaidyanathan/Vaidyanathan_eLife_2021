%% %%%%% FIGURE 6 SUPPLEMENT 2 %%%%%%

%% Data sets used in this figure
% Gq DREADD All
pathAll = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\DREADDs\Gq DREADDs\Gq-All'; %change this to your path on your computer

%Load Gq-DREADD All
load(strcat(pathAll,'\ephysData'))

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
  
%% Fig6FigSupp1 A

withBaseline = Gq(3:end); %these mice were injected immediately before recording and have an "elevated period"

split = 480; %in sec, time to split "elevated" vs "suppressed" in mice injected on wheel

%%% Figure Fig6FigSupp1A LEFT, sleep bout length

miceToUse = Gq;
sleepBoutLengthSal_1 = []; sleepBoutLengthSal_2 = [];
sleepBoutLengthCNO_1 = []; sleepBoutLengthCNO_2 = [];
for m = miceToUse
    if ismember(m, withBaseline) && ismember(m, lfpMice) 
        d = Sal_1{m};
        sleeps = allData(d).INT_lfp{2};
        sleeps1 = sleeps(sleeps(:,2) < split,:);
        sleeps2 = sleeps(sleeps(:,1) > split,:);
        hr1sLengths_pre = [sleeps1(:,2) - sleeps1(:,1)]';
        hr1sLengths_post = [sleeps2(:,2) - sleeps2(:,1)]';
        hr2sLengths = [allData(Sal_2{m}).INT_lfp{2}(:,2)-allData(Sal_2{m}).INT_lfp{2}(:,1)]';
        
        sleepBoutLengthSal_1 = [sleepBoutLengthSal_1 mean(hr1sLengths_pre)];
        sleepBoutLengthSal_2 = [sleepBoutLengthSal_2 mean([hr1sLengths_post hr2sLengths])];
               
        d = CNO_1{m};
        sleeps = allData(d).INT_lfp{2};
        sleeps1 = sleeps(sleeps(:,2) < split,:);
        sleeps2 = sleeps(sleeps(:,1) > split,:);
        hr1sLengths_pre = [sleeps1(:,2) - sleeps1(:,1)]';
        hr1sLengths_post = [sleeps2(:,2) - sleeps2(:,1)]';
        hr2sLengths = [allData(CNO_2{m}).INT_lfp{2}(:,2)-allData(CNO_2{m}).INT_lfp{2}(:,1)]';
        
        sleepBoutLengthCNO_1 = [sleepBoutLengthCNO_1 mean(hr1sLengths_pre)];
        sleepBoutLengthCNO_2 = [sleepBoutLengthCNO_2 mean([hr1sLengths_post hr2sLengths])];
         
    elseif ismember(m,lfpMice)       
       d = Sal_1{m};
        sleeps = allData(d).INT_lfp{2};
        hr1sLengths = [sleeps(:,2) - sleeps(:,1)]';
        hr2sLengths = [allData(Sal_2{m}).INT_lfp{2}(:,2)-allData(Sal_2{m}).INT_lfp{2}(:,1)]';
        sleepBoutLengthSal_2 = [sleepBoutLengthSal_2 mean([hr1sLengths hr2sLengths])];


        d = CNO_1{m};
        sleeps = allData(d).INT_lfp{2};
        hr1sLengths = [sleeps(:,2) - sleeps(:,1)]';
        hr2sLengths = [allData(CNO_2{m}).INT_lfp{2}(:,2)-allData(CNO_2{m}).INT_lfp{2}(:,1)]';
        sleepBoutLengthCNO_2 = [sleepBoutLengthCNO_2 mean([hr1sLengths hr2sLengths])];
    
    end    
end

%make NaNs = 0
sleepBoutLengthSal_1(isnan(sleepBoutLengthSal_1)) = 0;
sleepBoutLengthCNO_1(isnan(sleepBoutLengthCNO_1)) = 0;
sleepBoutLengthSal_2(isnan(sleepBoutLengthSal_2)) = 0;
sleepBoutLengthCNO_2(isnan(sleepBoutLengthCNO_2)) = 0;

figure;
clear s; k = 1;
for m = 1:length(sleepBoutLengthSal_2)
    s(k) = scatter([1 2], [sleepBoutLengthSal_2(m) sleepBoutLengthCNO_2(m)],80,'filled');
    hold on;
    errorbar([1 2], [sleepBoutLengthSal_2(m) sleepBoutLengthCNO_2(m)], [0 0],'k');
    hold on;
    xlim([.5 2.5]);
    k = k+1;
end
scatter([.75 2.25], [mean(sleepBoutLengthSal_2) mean(sleepBoutLengthCNO_2)],[],'k','filled');
hold on;
errorbar([.75 2.25], [mean(sleepBoutLengthSal_2) mean(sleepBoutLengthCNO_2)], [std(sleepBoutLengthSal_2) std(sleepBoutLengthCNO_2)],...
        'k','linestyle','none') 
[h,p] = ttest(sleepBoutLengthSal_2, sleepBoutLengthCNO_2);
    hold on;
text(1.2, max([sleepBoutLengthSal_2 sleepBoutLengthCNO_2]), strcat('p =',{' '},num2str(p)))  
xticks([1 2]); xticklabels({'saline', 'CNO'});
ylabel('mean NREM sleep bout duration'); %ylim([0 50])
title(strcat('"Low Ca2+" period: after',{' '}, num2str(split),{' '},'s'));
set(gca,'fontsize',15)
set(gcf,'renderer','painters')

%%% Figure Fig6FigSupp1A RIGHT, num of sleep bouts

miceToUse = Gq;
numSleepBoutsSal_1 = []; numSleepBoutsSal_2 = [];
numSleepBoutsCNO_1 = []; numSleepBoutsCNO_2 = [];
for m = miceToUse
    if ismember(m, withBaseline) && ismember(m, lfpMice) 
        d = Sal_1{m};
        sleeps = allData(d).INT_lfp{2};
        sleeps1 = sleeps(sleeps(:,2) < split,:);
        sleeps2 = sleeps(sleeps(:,1) > split,:);
        hr1sNum_pre = size(sleeps1,1);
        hr1sNum_post = size(sleeps2,1);
        hr2sNum = size(allData(Sal_2{m}).INT_lfp{2},1);
        
        numSleepBoutsSal_1 = [numSleepBoutsSal_1 hr1sNum_pre];
        numSleepBoutsSal_2 = [numSleepBoutsSal_2 sum([hr1sNum_post hr2sNum])];
               
        d = CNO_1{m};
        sleeps = allData(d).INT_lfp{2};
        sleeps1 = sleeps(sleeps(:,2) < split,:);
        sleeps2 = sleeps(sleeps(:,1) > split,:);
        hr1sNum_pre = size(sleeps1,1);
        hr1sNum_post = size(sleeps2,1);
        hr2sNum = size(allData(CNO_2{m}).INT_lfp{2},1);
        
        numSleepBoutsCNO_1 = [numSleepBoutsCNO_1 hr1sNum_pre];
        numSleepBoutsCNO_2 = [numSleepBoutsCNO_2 sum([hr1sNum_post hr2sNum])];
           
    elseif ismember(m,lfpMice)       
       d = Sal_1{m};
        sleeps = allData(d).INT_lfp{2};
        hr1sNum = size(sleeps,1);
        hr2sNum = size(allData(Sal_2{m}).INT_lfp{2},1);
        numSleepBoutsSal_2 = [numSleepBoutsSal_2 sum([hr1sNum hr2sNum])];


        d = CNO_1{m};
        sleeps = allData(d).INT_lfp{2};
        hr1sNum = size(sleeps,1);
        hr2sNum = size(allData(CNO_2{m}).INT_lfp{2},1);
        numSleepBoutsCNO_2 = [numSleepBoutsCNO_2 sum([hr1sNum hr2sNum])];
    
    end    
end

%make NaNs = 0
numSleepBoutsSal_1(isnan(numSleepBoutsSal_1)) = 0;
numSleepBoutsCNO_1(isnan(numSleepBoutsCNO_1)) = 0;
numSleepBoutsSal_2(isnan(numSleepBoutsSal_2)) = 0;
numSleepBoutsCNO_2(isnan(numSleepBoutsCNO_2)) = 0;    

figure;
clear s; k = 1;
for m = 1:length(numSleepBoutsSal_2)
    s(k) = scatter([1 2], [numSleepBoutsSal_2(m) numSleepBoutsCNO_2(m)],80,'filled');
    hold on;
    errorbar([1 2], [numSleepBoutsSal_2(m) numSleepBoutsCNO_2(m)], [0 0],'k');
    hold on;
    xlim([.5 2.5]);
    k = k+1;
end
scatter([.75 2.25], [mean(numSleepBoutsSal_2) mean(numSleepBoutsCNO_2)],[],'k','filled');
hold on;
errorbar([.75 2.25], [mean(numSleepBoutsSal_2) mean(numSleepBoutsCNO_2)], [std(numSleepBoutsSal_2) std(numSleepBoutsCNO_2)],...
        'k','linestyle','none') 
[h,p] = ttest(numSleepBoutsSal_2, numSleepBoutsCNO_2);
    hold on;
text(1.2, max([numSleepBoutsSal_2 numSleepBoutsCNO_2]), strcat('p =',{' '},num2str(p)))  
xticks([1 2]); xticklabels({'saline', 'CNO'});
ylabel('num of NREM sleep bouts'); %ylim([0 50])
title(strcat('"Low Ca2+" period: after',{' '}, num2str(split),{' '},'s'));
set(gca,'fontsize',15)
set(gcf,'renderer','painters')

%% Fig6FigSupp1 B

%%% Figure Fig6FigSupp1B LEFT, wake bout length

miceToUse = Gq;
wakeBoutLengthSal_1 = []; wakeBoutLengthSal_2 = [];
wakeBoutLengthCNO_1 = []; wakeBoutLengthCNO_2 = [];
for m = miceToUse
    if ismember(m, withBaseline) && ismember(m, lfpMice) 
        d = Sal_1{m};
        sleeps = allData(d).INT_lfp{1};
        sleeps1 = sleeps(sleeps(:,2) < split,:);
        sleeps2 = sleeps(sleeps(:,1) > split,:);
        hr1sLengths_pre = [sleeps1(:,2) - sleeps1(:,1)]';
        hr1sLengths_post = [sleeps2(:,2) - sleeps2(:,1)]';
        hr2sLengths = [allData(Sal_2{m}).INT_lfp{2}(:,2)-allData(Sal_2{m}).INT_lfp{2}(:,1)]';
        
        wakeBoutLengthSal_1 = [wakeBoutLengthSal_1 mean(hr1sLengths_pre)];
        wakeBoutLengthSal_2 = [wakeBoutLengthSal_2 mean([hr1sLengths_post hr2sLengths])];
               
        d = CNO_1{m};
        sleeps = allData(d).INT_lfp{1};
        sleeps1 = sleeps(sleeps(:,2) < split,:);
        sleeps2 = sleeps(sleeps(:,1) > split,:);
        hr1sLengths_pre = [sleeps1(:,2) - sleeps1(:,1)]';
        hr1sLengths_post = [sleeps2(:,2) - sleeps2(:,1)]';
        hr2sLengths = [allData(CNO_2{m}).INT_lfp{2}(:,2)-allData(CNO_2{m}).INT_lfp{2}(:,1)]';
        
        wakeBoutLengthCNO_1 = [wakeBoutLengthCNO_1 mean(hr1sLengths_pre)];
        wakeBoutLengthCNO_2 = [wakeBoutLengthCNO_2 mean([hr1sLengths_post hr2sLengths])];
         
    elseif ismember(m,lfpMice)       
       d = Sal_1{m};
        sleeps = allData(d).INT_lfp{1};
        hr1sLengths = [sleeps(:,2) - sleeps(:,1)]';
        hr2sLengths = [allData(Sal_2{m}).INT_lfp{2}(:,2)-allData(Sal_2{m}).INT_lfp{2}(:,1)]';
        wakeBoutLengthSal_2 = [wakeBoutLengthSal_2 mean([hr1sLengths hr2sLengths])];


        d = CNO_1{m};
        sleeps = allData(d).INT_lfp{1};
        hr1sLengths = [sleeps(:,2) - sleeps(:,1)]';
        hr2sLengths = [allData(CNO_2{m}).INT_lfp{2}(:,2)-allData(CNO_2{m}).INT_lfp{2}(:,1)]';
        wakeBoutLengthCNO_2 = [wakeBoutLengthCNO_2 mean([hr1sLengths hr2sLengths])];
    
    end    
end

%make NaNs = 0
wakeBoutLengthSal_1(isnan(wakeBoutLengthSal_1)) = 0;
wakeBoutLengthCNO_1(isnan(wakeBoutLengthCNO_1)) = 0;
wakeBoutLengthSal_2(isnan(wakeBoutLengthSal_2)) = 0;
wakeBoutLengthCNO_2(isnan(wakeBoutLengthCNO_2)) = 0;

figure;
clear s; k = 1;
for m = 1:length(wakeBoutLengthSal_2)
    s(k) = scatter([1 2], [wakeBoutLengthSal_2(m) wakeBoutLengthCNO_2(m)],80,'filled');
    hold on;
    errorbar([1 2], [wakeBoutLengthSal_2(m) wakeBoutLengthCNO_2(m)], [0 0],'k');
    hold on;
    xlim([.5 2.5]);
    k = k+1;
end
scatter([.75 2.25], [mean(wakeBoutLengthSal_2) mean(wakeBoutLengthCNO_2)],[],'k','filled');
hold on;
errorbar([.75 2.25], [mean(wakeBoutLengthSal_2) mean(wakeBoutLengthCNO_2)], [std(wakeBoutLengthSal_2) std(wakeBoutLengthCNO_2)],...
        'k','linestyle','none') 
[h,p] = ttest(wakeBoutLengthSal_2, wakeBoutLengthCNO_2);
    hold on;
text(1.2, max([wakeBoutLengthSal_2 wakeBoutLengthCNO_2]), strcat('p =',{' '},num2str(p)))  
xticks([1 2]); xticklabels({'saline', 'CNO'});
ylabel('mean wake bout duration'); ylim([0 500])
title(strcat('"Low Ca2+" period: after',{' '}, num2str(split),{' '},'s'));
set(gca,'fontsize',15)
set(gcf,'renderer','painters')

%%% Figure Fig6FigSupp1B RIGHT, num of wake bouts

miceToUse = Gq;
numWakeBoutsSal_1 = []; numWakeBoutsSal_2 = [];
numWakeBoutsCNO_1 = []; numWakeBoutsCNO_2 = [];
for m = miceToUse
    if ismember(m, withBaseline) && ismember(m, lfpMice) 
        d = Sal_1{m};
        sleeps = allData(d).INT_lfp{1};
        sleeps1 = sleeps(sleeps(:,2) < split,:);
        sleeps2 = sleeps(sleeps(:,1) > split,:);
        hr1sNum_pre = size(sleeps1,1);
        hr1sNum_post = size(sleeps2,1);
        hr2sNum = size(allData(Sal_2{m}).INT_lfp{1},1);
        
        numWakeBoutsSal_1 = [numWakeBoutsSal_1 hr1sNum_pre];
        numWakeBoutsSal_2 = [numWakeBoutsSal_2 sum([hr1sNum_post hr2sNum])];
               
        d = CNO_1{m};
        sleeps = allData(d).INT_lfp{1};
        sleeps1 = sleeps(sleeps(:,2) < split,:);
        sleeps2 = sleeps(sleeps(:,1) > split,:);
        hr1sNum_pre = size(sleeps1,1);
        hr1sNum_post = size(sleeps2,1);
        hr2sNum = size(allData(CNO_2{m}).INT_lfp{1},1);
        
        numWakeBoutsCNO_1 = [numWakeBoutsCNO_1 hr1sNum_pre];
        numWakeBoutsCNO_2 = [numWakeBoutsCNO_2 sum([hr1sNum_post hr2sNum])];
           
    elseif ismember(m,lfpMice)       
       d = Sal_1{m};
        sleeps = allData(d).INT_lfp{1};
        hr1sNum = size(sleeps,1);
        hr2sNum = size(allData(Sal_2{m}).INT_lfp{1},1);
        numWakeBoutsSal_2 = [numWakeBoutsSal_2 sum([hr1sNum hr2sNum])];


        d = CNO_1{m};
        sleeps = allData(d).INT_lfp{1};
        hr1sNum = size(sleeps,1);
        hr2sNum = size(allData(CNO_2{m}).INT_lfp{1},1);
        numWakeBoutsCNO_2 = [numWakeBoutsCNO_2 sum([hr1sNum hr2sNum])];
    
    end    
end

%make NaNs = 0
numWakeBoutsSal_1(isnan(numWakeBoutsSal_1)) = 0;
numWakeBoutsCNO_1(isnan(numWakeBoutsCNO_1)) = 0;
numWakeBoutsSal_2(isnan(numWakeBoutsSal_2)) = 0;
numWakeBoutsCNO_2(isnan(numWakeBoutsCNO_2)) = 0;

figure;
clear s; k = 1;
for m = 1:length(numWakeBoutsSal_2)
    s(k) = scatter([1 2], [numWakeBoutsSal_2(m) numWakeBoutsCNO_2(m)],80,'filled');
    hold on;
    errorbar([1 2], [numWakeBoutsSal_2(m) numWakeBoutsCNO_2(m)], [0 0],'k');
    hold on;
    xlim([.5 2.5]);
    k = k+1;
end
scatter([.75 2.25], [mean(numWakeBoutsSal_2) mean(numWakeBoutsCNO_2)],[],'k','filled');
hold on;
errorbar([.75 2.25], [mean(numWakeBoutsSal_2) mean(numWakeBoutsCNO_2)], [std(numWakeBoutsSal_2) std(numWakeBoutsCNO_2)],...
        'k','linestyle','none') 
[h,p] = ttest(numWakeBoutsSal_2, numWakeBoutsCNO_2);
    hold on;
text(1.2, max([numWakeBoutsSal_2 numWakeBoutsCNO_2]), strcat('p =',{' '},num2str(p)))  
xticks([1 2]); xticklabels({'saline', 'CNO'});
ylabel('num of wake bouts'); ylim([-10 100])
title(strcat('"Low Ca2+" period: after',{' '}, num2str(split),{' '},'s'));
set(gca,'fontsize',15)
set(gcf,'renderer','painters') 
