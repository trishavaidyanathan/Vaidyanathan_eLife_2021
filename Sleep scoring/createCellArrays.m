function [uniqueMice, mice, CNO, Sal, CNO_1, CNO_2, CNO_3, Sal_1, Sal_2, Sal_3, hasbaseline,...
    lfpMice, eegMice, lfpFiles, eegFiles] = createCellArrays(allData)
%% Creates cell arrays corresponding to files for each mouse

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

lfpMice = []; %corresponds to uniqueMice, which mice have lfp rec
eegMice = []; %corresponds to uniqueMice, which mice have eeg rec
lfpFiles = [];
eegFiles = [];
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

if contains(allData(9).dreadd,'Gi') && ~contains(allData(9).dreadd,'IP3R2KO')
    Sal_2{9} = Sal_1{9}(1);
    Sal_1{9} = Sal_1{9}(2); %%% TEMPORARY for Gi files
    disp('changed sal hour 2 file')
end

hasbaseline = [];
for m = 1:length(uniqueMice)
    currSal = Sal_1{m}; currCNO = CNO_1{m};
    if ~isempty(allData(currSal).dataBaseline) && ~isempty(allData(currCNO).dataBaseline)
        hasbaseline = [hasbaseline m];
    end
end
end