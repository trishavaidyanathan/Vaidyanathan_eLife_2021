%% %%%%% FIGURE 6 FIGURE SUPPLEMENT 1 %%%%%%

%% Data sets used in this figure
% Immuno colocalization data
immunopath = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\Immunostaining\Gq\cellCounts_Gq'; %change to the file path on your computer for the Gi-immuno XLS file
% WEKA Gi Data
wekapath = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\Immunostaining\Gq'; %change this to your path on your computer
% Gq DREADD CNO doses
pathDoses = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\DREADDs\Gq DREADDs\Compare CNO doses'; %change this to your path on your computer
% ROI data (soma + processes)
pathROI = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\DREADDs\Gq DREADDs\Soma-Proccesses ROIs'; %change this to your path on your computer
% Wake cocktail data
pathCocktail = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\Wake Cocktail'; %change this to your path on your computer

%Load ROI data
load(strcat(pathROI,'\soma-processes_ROIs'));

%Load Weka Data:
load(strcat(wekapath,'\','wekaGq'));

%Load Gq doses data 
load(strcat(pathDoses,'\ephysData'))
load(strcat(pathDoses,'\aquaData'))

%Load Wake cocktail data
load(strcat(pathCocktail,'\wakeCocktail_AQuA'));

%% Pre-processing FOR Gq CNO DOSE DATA (Assumes data is in "allData" and "allRes") 
    %% Create uniqueMice cell-array
        sf=1000;

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

            %Sort saline condition into hours (regardless of day)
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
        lfpFiles = []; %corresponds to allData, which files have lfp rec
        eegMice = []; %corresponds to uniqueMice, which mice have eeg rec
        eegFiles = []; %corresponds to allData, which files have eeg rec
        for m = 1:length(uniqueMice)
            if ~isempty(allData(Sal_1{m}(1)).lfp) && ~isempty(allData(CNO_1{m}(1)).lfp)
                lfpMice = [lfpMice m];
        %         lfpFiles = [lfpFiles mice{m}'];
            end
            if ~isempty(allData(Sal_1{m}(1)).eeg) && ~isempty(allData(CNO_1{m}(1)).eeg)
                eegMice = [eegMice m];
        %         eegFiles = [eegFiles mice{m}'];
            end
            if ~isempty(allData(Sal_1{m}(1)).lfp)
                lfpFiles = [lfpFiles [Sal_1{m}]];
            end
            if ~isempty(allData(Sal_1{m}(1)).eeg)
                eegFiles = [eegFiles [Sal_1{m}]];
            end
            if ~isempty(allData(CNO_1{m}(1)).lfp)
                lfpFiles = [lfpFiles [CNO_1{m}]];
            end
            if ~isempty(allData(CNO_1{m}(1)).eeg)
                eegFiles = [eegFiles [CNO_1{m}]];
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

%% Fig6FigSupp1 C
%%% USES COLOCALIZATION DATA, NOTHING TO LOAD

numSheets = 12; %corresponding to number of mice (1 mouse per sheet)

for i = 1:numSheets 
    [data,txt] = xlsread(immunopath,i);
    colocData(i).name = txt{1};
    colocData(i).data = data;
end

%%% REMOVE following mice due to staining issues 
colocData(1) = []; %remove m116
colocData(1) = []; %remove m117

percentColoc1 = [];
percentColoc2 = [];
for i = 1:length(colocData)
    percentColoc1 = [percentColoc1 nanmean(colocData(i).data(:,2)./colocData(i).data(:,1))*100];
    summed = nansum(colocData(i).data);
    percentColoc2 = [percentColoc2 (summed(2)/summed(1))*100];
end

figure;
bar(mean(percentColoc1)); 
hold on;
s1 = scatter(ones(1,length(percentColoc1)), percentColoc1,'filled', 'jitter','on','jitterAmount',.15);
s1.MarkerFaceAlpha = .5;
ylabel('% NeuN + mCherry+')
title('Gq colocalization neurons + mCherry');
xlim([-1 3]); ylim([0 100])
set(gca,'fontsize',15)
text(1,50,strcat('mean =',{' '},num2str(mean(percentColoc1))))
text(1,45,strcat('SD =',{' '},num2str(std(percentColoc1))))
text(1,40,strcat('SEM =',{' '},num2str(sem(percentColoc1))))

%% Fig6FigSupp1 E
%%% LOAD WEKA DATA

%for allen brain atlas slice numbers
sliceDiff = .1; %in mm, from Allen Atlas, each slice num is 100 um apart
bregmaDist = fliplr([0:132]*sliceDiff); %bregma distances from allen atlas
bregmaDist(81) = [];
bregmaDist(60:end) = bregmaDist(60:end) + -7.855;
bregmaDist(1:59) = bregmaDist(1:59) + -7.88; 

binwidth = 4; %bin images in 4 slice increments 
allBins = cellfun(@(x) x', {wekaData(:).bins},'UniformOutput',0);
allBins = cell2mat(allBins);
breadth = max(allBins)-min(allBins);
numbins = floor(breadth/binwidth);
binSt = (1:numbins)*binwidth + min(allBins);
binSt = [min(allBins) binSt(1:end-1)];

idx = repmat(0:binwidth-1,numbins,1) + binSt';

ipsiPixAll = cellfun(@(x) x', {wekaData(:).IPSIpix},'UniformOutput',0);
ipsiPixAll = cell2mat(ipsiPixAll)';
contraPixAll = cellfun(@(x) x', {wekaData(:).CONTRApix},'UniformOutput',0);
contraPixAll = cell2mat(contraPixAll)';
ipsiVals = {}; contraVals = {};
for b = 1:length(binSt)
    currBins = find(ismember(allBins,idx(b,:)));
    ipsiVals{b} = ipsiPixAll(currBins,2)./ipsiPixAll(currBins,1);
    contraVals{b} = contraPixAll(currBins,2)./ipsiPixAll(currBins,1);
end

bin_mm = bregmaDist(idx);
binNames = mean(bin_mm,2); %bin center in mm from bregma

refScrew = 84; %From excel sheet, found average allen slice # of skull screw from each slice image
%Gq: 84, Gi: 89;
refScrewMM = bregmaDist(refScrew);
FCscrew = 2.7; 

%%% Plot box plot

figure;
subplot(2,1,1)
ipsiL = cellfun(@length, ipsiVals);
g = []; for b=1:length(binSt); g = [g b*ones(1,ipsiL(b))]; end 
ipsi = cellfun(@(x) x', ipsiVals,'UniformOutput',0);
boxplot([ipsi{:}],g,'positions',binNames,'labels',binNames)
hold on; scatter(refScrewMM, .3, 60, 'r','*'); text(refScrewMM,.3,'LFP');
set(gca,'Xdir','reverse')
hold on;
for b = 1:length(binSt)
    vals = ipsiVals{b};
    vals = vals(~isnan(vals));
    s = scatter(ones(1,length(vals))*binNames(b), vals, 'r', 'filled', 'jitter','on','jitterAmount',.15);
    hold on;
end
ylim([-.05 .5]);
ylabel('spread of expression'); xlabel('mm from bregma');
title('ipsilateral to injection site');
set(gca,'fontsize',15)
set(gcf,'renderer','painters');
xlim([-4.8 2.8])

%Bar plot for CONTRA
subplot(2,1,2)
contraL = cellfun(@length, contraVals);
g = []; for b=1:length(binSt); g = [g b*ones(1,contraL(b))]; end 
contra = cellfun(@(x) x', contraVals,'UniformOutput',0);
boxplot([contra{:}],g,'positions',binNames,'labels',binNames)
set(gca,'Xdir','reverse')
hold on; scatter(FCscrew, .3, 60, 'r','*'); text(FCscrew,.3,'FC-EEG');
hold on;
for b = 1:length(binSt)
    vals = contraVals{b};
    vals = vals(~isnan(vals));
    s = scatter(ones(1,length(vals))*binNames(b), vals, 'r', 'filled', 'jitter','on','jitterAmount',.15);
    s.MarkerFaceAlpha = .5;
    hold on;
end
ylim([-.05 .5]);xlim([-4.8 2.8])
ylabel('spread of expression'); xlabel('mm from bregma');
title('ipsilateral to injection site');
set(gca,'fontsize',15)
set(gcf,'renderer','painters');

%% Fig6FigSupp1 F 
%%% LOAD GQ DOSE DATA

injHeadFix = [5:18 27:38]; %need to split these files at split num (480)

%PERCENT DECREASE
dose01 = find(~cellfun(@isempty, strfind({allData(:).dose}, '0.1')));
dose05 = find(~cellfun(@isempty, strfind({allData(:).dose}, '0.5')));
dose1 = find(~cellfun(@isempty, strfind({allData(:).dose}, '1 mg/kg')) & ...
    cellfun(@isempty, strfind({allData(:).dose}, '0')));

dose01(1:2) = []; %remove extra m124 files


evtSal_01 = []; evtCNO_01 = [];
evtSal_05 = []; evtCNO_05 = [];
evtSal_1 = []; evtCNO_1 = [];
for m = 1:length(uniqueMice)
    currDose01 = intersect(mice{m},dose01);
    currDose05 = intersect(mice{m},dose05);
    currDose1 = intersect(mice{m},dose1);
    
    sal01 = intersect(currDose01,Sal{m});
    cno01 = intersect(currDose01,CNO{m});
    
    sal05 = intersect(currDose05,Sal{m});
    cno05 = intersect(currDose05,CNO{m});
    
    sal1 = intersect(currDose1,Sal{m});
    cno1 = intersect(currDose1,CNO{m});
    
    % DOSE 0.1 mg/kg
    %Saline
    d = sal01; if ~isempty(d) && ismember(d(1),injHeadFix) && ismember(d(2),injHeadFix)
        cumEvts_sal = cumsum([allRes(d(1)).numEvents allRes(d(2)).numEvents]);
        toFill = 7200-length(cumEvts_sal);
        cumEvts_sal = [cumEvts_sal NaN(1, toFill)];  
        maxSal = max(cumEvts_sal); 
        %CNO
        d = cno01; 
        cumEvts_cno = cumsum([allRes(d(1)).numEvents allRes(d(2)).numEvents]);
        toFill = 7200-length(cumEvts_cno);
        cumEvts_cno = [cumEvts_cno NaN(1, toFill)]; 

        evtSal_01 = [evtSal_01 ; (cumEvts_sal/maxSal)*100];
        evtCNO_01 = [evtCNO_01 ; (cumEvts_cno/maxSal)*100];
    end 
    
    % DOSE 0.5 mg/kg
    %Saline
    d = sal05; if ~isempty(d) && ismember(d(1),injHeadFix) && ismember(d(2),injHeadFix)
        cumEvts_sal = cumsum([allRes(d(1)).numEvents allRes(d(2)).numEvents]);
        toFill = 7200-length(cumEvts_sal);
        cumEvts_sal = [cumEvts_sal NaN(1, toFill)]; 
        maxSal = max(cumEvts_sal);  
        %CNO
        d = cno05; 
        cumEvts_cno = cumsum([allRes(d(1)).numEvents allRes(d(2)).numEvents]);
        toFill = 7200-length(cumEvts_cno);
        cumEvts_cno = [cumEvts_cno NaN(1, toFill)]; 

        evtSal_05 = [evtSal_05 ; (cumEvts_sal/maxSal)*100];
        evtCNO_05 = [evtCNO_05 ; (cumEvts_cno/maxSal)*100];
    end
    
    % DOSE 1 mg/kg
    %Saline
    d = sal1; if ~isempty(d) && ismember(d(1),injHeadFix) && ismember(d(2),injHeadFix)
        cumEvts_sal = cumsum([allRes(d(1)).numEvents allRes(d(2)).numEvents]);
        toFill = 7200-length(cumEvts_sal);
        cumEvts_sal = [cumEvts_sal NaN(1, toFill)]; 
        maxSal = max(cumEvts_sal); 
        %CNO
        d = cno1;
        cumEvts_cno = cumsum([allRes(d(1)).numEvents allRes(d(2)).numEvents]);
        toFill = 7200-length(cumEvts_cno);
        cumEvts_cno = [cumEvts_cno NaN(1, toFill)]; 

        evtSal_1 = [evtSal_1 ; (cumEvts_sal/maxSal)*100];
        evtCNO_1 = [evtCNO_1 ; (cumEvts_cno/maxSal)*100];
    end
end

figure;
    meanSal = nanmean(evtSal_01,1);
    meanCNO = nanmean(evtCNO_01,1);
    semSal = nanstd(evtSal_01) / sqrt(size(evtSal_01,1));
    semCNO = nanstd(evtCNO_01) / sqrt(size(evtCNO_01,1));
    shadedErrorBar([1:length(meanSal)]/60, meanSal, semSal, ...
        'lineprops','-k','transparent',1)
    hold on;
    shadedErrorBar([1:length(meanCNO)]/60, meanCNO, semCNO, ...
        'lineprops','-r','transparent',1)    
    ylabel('cumulative events (% saline total)');
    xlabel('time (min)');
    title('0.1 mg/kg CNO');
    xlim([0 120])
    hold on; vline(19,'--k');
    set(gca,'fontsize',15)

%% Fig6FigSupp1 G
%%% LOAD GQ DOSE DATA

%%% Uses same variables (evtSal_1, evtCNO_1) as section above
    meanSal = nanmean(evtSal_1,1);
    meanCNO = nanmean(evtCNO_1,1);
    semSal = nanstd(evtSal_1) / sqrt(size(evtSal_1,1));
    semCNO = nanstd(evtCNO_1) / sqrt(size(evtCNO_1,1));
figure;
    shadedErrorBar([1:length(meanSal)]/60, meanSal, semSal, ...
        'lineprops','-k','transparent',1)
    hold on;
    shadedErrorBar([1:length(meanCNO)]/60, meanCNO, semCNO, ...
        'lineprops','-r','transparent',1)    
    ylabel('cumulative events (% saline total)');
    xlabel('time (min)');
    title('1 mg/kg CNO');
    xlim([0 120])
    hold on; vline(8,'--k');
    set(gca,'fontsize',15)
    set(gcf,'renderer','painters')

%% Fig6FigSupp1 H
%%% LOAD GQ DOSE DATA

percentDecr_01 = {};
percentDecr_05 = {};
percentDecr_1 = {};
hourCountSal = zeros(3,length(uniqueMice)); hourCountCNO = zeros(3,length(uniqueMice)); %for "n" for paper, rows correspond to 0.1, 0.5, 1 mg/kg 
for m = 1:length(uniqueMice)
    currDose01 = intersect(mice{m},dose01);
    currDose05 = intersect(mice{m},dose05);
    currDose1 = intersect(mice{m},dose1);
    
    sal01 = intersect(currDose01,Sal{m});
    cno01 = intersect(currDose01,CNO{m});
    
    sal05 = intersect(currDose05,Sal{m});
    cno05 = intersect(currDose05,CNO{m});
    
    sal1 = intersect(currDose1,Sal{m});
    cno1 = intersect(currDose1,CNO{m});
    
    hourCountSal(1,m) = length(sal01);
    hourCountSal(2,m) = length(sal05);
    hourCountSal(3,m) = length(sal1);
    hourCountCNO(1,m) = length(cno01);
    hourCountCNO(2,m) = length(cno05);
    hourCountCNO(3,m) = length(cno1);
    
    split = 1200;
    temp = {allRes([sal01]').numEvents};
    if ismember(sal01, injHeadFix); temp = cellfun(@(x) x(split+1:end), temp,'UniformOutput',0); end 
    if length(temp) == 1 || length(temp) > 2; disp('not 2 hrs of data in this condition'); end
    sal01_ER = sum(cellfun(@sum, temp,'UniformOutput',1));
    
    temp = {allRes([cno01]').numEvents};
    if ismember(cno01, injHeadFix); temp = cellfun(@(x) x(split+1:end), temp,'UniformOutput',0); end 
    if length(temp) == 1 || length(temp) > 2; disp('not 2 hrs of data in this condition'); end
    cno01_ER = sum(cellfun(@sum, temp,'UniformOutput',1));
    
    temp = {allRes([sal05]').numEvents};
    if ismember(sal05, injHeadFix); temp = cellfun(@(x) x(split+1:end), temp,'UniformOutput',0); end 
    if length(temp) == 1 || length(temp) > 2; disp('not 2 hrs of data in this condition'); end
    sal05_ER = sum(cellfun(@sum, temp,'UniformOutput',1));
    
    temp = {allRes([cno05]').numEvents};
    if ismember(cno05, injHeadFix); temp = cellfun(@(x) x(split+1:end), temp,'UniformOutput',0); end 
    if length(temp) == 1 || length(temp) > 2; disp('not 2 hrs of data in this condition'); end
    cno05_ER = sum(cellfun(@sum, temp,'UniformOutput',1));
    
    split = 1200;
    temp = {allRes([sal1]').numEvents};
    if ismember(sal1, injHeadFix); temp = cellfun(@(x) x(split+1:end), temp,'UniformOutput',0); end 
    if length(temp) == 1 || length(temp) > 2; disp('not 2 hrs of data in this condition'); end
    sal1_ER = sum(cellfun(@sum, temp,'UniformOutput',1));
    
    temp = {allRes([cno1]').numEvents};
    if ismember(cno1, injHeadFix); temp = cellfun(@(x) x(split+1:end), temp,'UniformOutput',0); end 
    if length(temp) == 1 || length(temp) > 2; disp('not 2 hrs of data in this condition'); end
    cno1_ER = sum(cellfun(@sum, temp,'UniformOutput',1)); 
    
    percentDecr_01{m} = ((cno01_ER-sal01_ER)/sal01_ER)*100;
    percentDecr_05{m} = ((cno05_ER-sal05_ER)/sal05_ER)*100;
    percentDecr_1{m} = ((cno1_ER-sal1_ER)/sal1_ER)*100;
end


figure;
scatter([1 2 3], [nanmean([percentDecr_01{:}]), nanmean([percentDecr_05{:}]), nanmean([percentDecr_1{:}])],'k','filled')
hold on;
errorbar([1 2 3], [nanmean([percentDecr_01{:}]), nanmean([percentDecr_05{:}]), nanmean([percentDecr_1{:}])],...
    [sem([percentDecr_01{:}]), sem([percentDecr_05{:}]), sem([percentDecr_1{:}])],'linestyle','none')
xlim([0 4])
xticks([1 2 3]); xticklabels({'0.1 mg/kg','0.5 mg/kg','1 mg/kg'})
ylim([-100 10])
hline(0,'--k')
ylabel('% change in event rate')
title('split at 1200s')
set(gca,'fontsize',15);

%% Fig6FigSupp1 I 
%%% LOAD ROI DATA

% CALCULATE DFOF

win = 100; %frames;
for d = 1:length(allData)
    baselineP = mean(allData(d).tracesP(:,1:win),2); %baseline for each ROI
    dfofP = (allData(d).tracesP - baselineP) ./ baselineP;
    allData(d).dfofPall = dfofP;
    allData(d).baselineP = baselineP;
    
    baselineS = mean(allData(d).tracesS(:,1:win),2); %baseline for each ROI
    dfofS = (allData(d).tracesS - baselineS) ./ baselineS;
    allData(d).dfofS = dfofS;
    allData(d).baselineS = baselineS;
end

% REMOVE P-ROIs WITHOUT SIGNAL

allthresh = [];
for d = 1:length(allData)
    allMags = []; allLocs = [];
    amp = max(allData(d).dfofPall') - mean(allData(d).dfofPall(:,1:win),2)'; 
    thresh = mean(amp);
    allthresh = [allthresh thresh];
        
    withSignal = [];
    for r = 1:size(allData(d).dfofPall,1) %for each ROI
        curr = allData(d).dfofPall(r,:);
        [l,m] = peakfinder(smooth(curr,15),thresh,thresh,1,0,0); %sel, thresh, extrema, endpoints, interpolate
        [mag,i] = max(m);
        loc = l(i);
        if ~isempty(mag)
            allMags = [allMags mag];
            allLocs = [allLocs loc];
            withSignal = [withSignal r]; %add this ROI to list
        end
    end
    allData(d).signalMagsP = allMags;
    allData(d).signalLocsP = allLocs; 
    allData(d).dfofP = allData(d).dfofPall(withSignal,:); %just select out ROIs with signal
end

% FIND FRAMETIMES

for d = 1:length(allData)
    allData(d).frameperiod = 0.6034953869;
    allData(d).frameTimes = [1:size(allData(d).tracesS,2)]*0.6034953869;
end

% PLOT

allSomas = [];
allPs = [];
for d = 1:length(allData)
    allSomas = [allSomas ; allData(d).dfofS ./ max(allData(d).dfofS')'];
    allPs = [allPs ; allData(d).dfofP ./ max(allData(d).dfofP')'];
end
somaMean = mean(allSomas);
somaSD = std(allSomas) ./ sqrt(size(allSomas,1));
pMean = mean(allPs);
pSD = std(allPs) ./ sqrt(size(allPs,1));

figure;
s1=shadedErrorBar(allData(d).frameTimes/60, somaMean, somaSD,'lineprops',{'-k','markerfacecolor',[0 0.45 0.75]});
hold on;
s2=shadedErrorBar(allData(d).frameTimes/60, pMean, pSD,'lineprops',{'-r','markerfacecolor',[0.85 0.33 0.099]});
legend([s1.mainLine s2.mainLine],{'somas','processes'})
xlim([0 60])
xlabel('time (min)');
ylabel('dfof (normalized to peak)');
title('Gq CNO, n = 3 mice')
set(gca,'fontsize',15)
set(gcf,'renderer','painters')

%% Fig6FigSupp1 J
%%% LOAD WAKE COCKTAIL DATA

for d = 1:length(allData)
    if ~isempty(allData(d).opts)
        allData(d).frameperiod = allData(d).opts.frameRate;
        allData(d).numFrames = allData(d).opts.sz(3);
        allData(d).peakTimesF = allData(d).fts.loc.t0; %in frames (onsets)
        allData(d).peakTimes = allData(d).peakTimesF*allData(d).frameperiod; %in sec
        [allData(d).numEvents, allData(d).timebins] = hist(allData(d).peakTimes, 1:(allData(d).numFrames * allData(d).frameperiod)); %
        allData(d).eventRate = allData(d).numEvents/(allData(d).timebins(2)-allData(d).timebins(1)); %events per second
    end
end

%find rows of allData corresponding to x1 dose of cocktail and x2 dose of
%cocktail
first = cellfun(@isempty,strfind([allData(:).name],'x2'));
second = ~first;

st = 80; %start st sec before 0 (cocktail add)
cumEvents1 = NaN(length(find(first)),length(-st:300));
k = 1;
for d = find(first)
    cocktailAdd = round(allData(d).cocktail1F * allData(d).frameperiod)+60;
    events = [NaN(1,st-cocktailAdd+1) cumsum(allData(d).numEvents)];
    cumEvents1(k,1:length(events)) = events;
    k = k+1;
end

st = 80; %start st sec before 0 (cocktail add)
cumEvents2 = NaN(length(find(second)),length(-st:300));
k = 1;
for d = find(second)
    cocktailAdd = round(allData(d).cocktail2F * allData(d).frameperiod)+60;
    events = [NaN(1,st-cocktailAdd+1) cumsum(allData(d).numEvents)];
    cumEvents2(k,1:length(events)) = events;
    k = k+1;
end
cumEvents2 = cumEvents2 + max(cumEvents1')';

norm1 = cumEvents1 ./ max(cumEvents2')';
norm2 = cumEvents2 ./ max(cumEvents2')';
minIdx1 = []; maxIdx1 = []; minIdx2 = []; maxIdx2 = [];
for i = 1:size(norm1,1); 
       minIdx1 = [minIdx1 min(find(~isnan(norm1(i,:))))]; 
       minIdx2 = [minIdx2 min(find(~isnan(norm2(i,:))))]; 
       maxIdx1 = [maxIdx1 max(find(~isnan(norm1(i,:))))]; 
       maxIdx2 = [maxIdx2 max(find(~isnan(norm2(i,:))))];
end
minIdx1 = max(minIdx1); minIdx2 = max(minIdx2); 
maxIdx1 = min(maxIdx1); maxIdx2 = min(maxIdx2);

figure; 
tIdx = [-st:300];
shadedErrorBar(tIdx(minIdx1:maxIdx1)/60, nanmean(norm1(:,minIdx1:maxIdx1)),nanstd(norm1(:,minIdx1:maxIdx1))/sqrt(size(norm1,1)))    
hold on;
vline(0,'--k'); text(0, .8, 'wake cocktail');
hold on;
tIdx = [-st:300]+(st+301);
shadedErrorBar(tIdx(minIdx2:maxIdx2)/60, nanmean(norm2(:,minIdx2:maxIdx2)),nanstd(norm2(:,minIdx2:maxIdx2))/sqrt(size(norm2,1)))    
hold on;
vline((st+301)/60,'--k'); text((st+301)/60, .8, 'wake cocktail x2');
ylabel('cumulative events (norm max)'); xlabel('time (min)');
set(gca,'fontsize',15)
set(gcf,'renderer','painters')
