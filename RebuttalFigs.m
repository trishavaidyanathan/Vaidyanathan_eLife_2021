%%%
%%%
%%%
% For rebuttal
%%%
%%%
%%%
%% Correlation bewteen EEG and LFP for SWA

% Using following mat file for saline:
load('E:\Analysis\MAT files\MAT Files for Paper 2020\DREADDs\All DREADDs Saline\dataAstroCa_Gq-and-Gi_SalineOnly.mat')

% Using following mat file for CNO
load('E:\Analysis\MAT files\MAT Files for Paper 2020\DREADDs\Gi DREADDs\dataDREADDs_Gi-allHours-wBaseline.mat')

freq = [0.5 4]; %use this frequency range for power
% not using specSelect because that's smoothed

allLFP = []; allEEG = []; %collect the average power in 'freq' for each file all in one vector
for d = 1:length(allData)
    if ~isempty(allData(d).lfp) && ~isempty(allData(d).eeg)
        [~,st] = min(abs(allData(d).f_FFT_lfp - freq(1)));
        [~,fn] = min(abs(allData(d).f_FFT_lfp - freq(2)));
        currLFP = zscore(mean(log(allData(d).FFT_lfp(:,st:fn)),2));
        
        allLFP = [allLFP ; currLFP];
        
        [~,st] = min(abs(allData(d).f_FFT_eeg - freq(1)));
        [~,fn] = min(abs(allData(d).f_FFT_eeg - freq(2)));
        currEEG = zscore(mean(log(allData(d).FFT_eeg(:,st:fn)),2));
        
        allEEG = [allEEG ; currEEG];
    end
end

[r,p] = corrcoef(allLFP,allEEG);

figure; scatter(allLFP,allEEG,'filled','MarkerFaceAlpha', .2)
text(-2, 6,strcat('r =',{' '},num2str(r(2))),'FontSize',15)
text(-2, 5.5, strcat('p =',{' '},num2str(p(2))),'FontSize',15)
xlabel('LFP SWA'); ylabel('EEG SWA');
title('Gi CNO')
set(gca,'fontsize',15)

%% Find time when sleep scoring agrees between LFP and EEG

percentAgree = [];
for d = 1:length(allData)
    if ~isempty(allData(d).lfp) && ~isempty(allData(d).eeg)
        lfpIDX = allData(d).IDX_lfp;
        eegIDX = allData(d).IDX_eeg;
        
        agree = (length(find(lfpIDX == eegIDX))/length(lfpIDX))*100;
        percentAgree = [percentAgree agree];
    end
end
mean(percentAgree)

%% Analyze dfof in soma/proccess ROIs with Gq 

% LOAD DATA

outputPath = 'E:\Analysis\Fiji Output Soma-Processes';

folders = dir(outputPath); 
folders = {folders(:).name}';
folders = folders(3:end);

filePaths = cellfun(@(x)[outputPath '\' x], folders,'uniformoutput',false);

d = 1;
for i = 1:length(filePaths)
    files = dir(filePaths{i});
    allData(d).mouseID = filePaths{i}(end-3:end);
    
    fileNames = lower({files(:).name});
    somas = ~cellfun(@isempty, strfind(fileNames,'soma')) & ~cellfun(@isempty, strfind(fileNames,'traces'));
    processes = ~cellfun(@isempty, strfind(fileNames,'process')) & ~cellfun(@isempty, strfind(fileNames,'traces'));
    
    somaPath = strcat(filePaths{i},'\',fileNames{somas});
    processPath = strcat(filePaths{i},'\',fileNames{processes});
    
    tracesS = csvread(somaPath,1,1);
    allData(d).tracesS = tracesS';
    tracesP = load(processPath);
    allData(d).tracesP = tracesP.traces;
    
    d = d + 1;
end
    

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

figOn = 0;
allthresh = [];
for d = 1:length(allData)
    allMags = []; allLocs = [];
    amp = max(allData(d).dfofPall') - mean(allData(d).dfofPall(:,1:win),2)'; 
    thresh = mean(amp);
    allthresh = [allthresh thresh];
    
    figure;
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

            if figOn
                plot(allData(d).dfofPall(r,:)+length(allMags))
                hold on;
                scatter(loc, mag + length(allMags), 'filled')
                hold on;
            end
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

%% Compare fluorescent intensiy for Gq vs Gi (63x colocalization data)
clear
GqPath = 'E:\Staining Results\DREAAD (Gq + LFP) Project\FluoIntensity_Colocal';
GiPath = 'E:\Staining Results\DREADD Gi\FluoIntensity_Colocal';

GqFiles = dir(GqPath);
GqFiles = GqFiles(3:end);
GiFiles = dir(GiPath);
GiFiles = GiFiles(3:end);

GqData_perMouse = {};
GqData_all = [];
for f = 1:length(GqFiles)
    GqData_perMouse{f} = csvread(strcat(GqPath,'\',GqFiles(f).name),1,1);
    GqData_all = [GqData_all ; GqData_perMouse{f}];
end

GiData_perMouse = {};
GiData_all = [];
for f = 1:length(GiFiles)
    GiData_perMouse{f} = csvread(strcat(GiPath,'\',GiFiles(f).name),1,1);
    GiData_all = [GiData_all ; GiData_perMouse{f}];
end

Gq = cellfun(@mean, GqData_perMouse);
Gi = cellfun(@mean, GiData_perMouse);

figure;
boxplot([Gi Gq],[repmat(1,1,length(Gi)),repmat(2,1,length(Gq))],'Labels',{'Gi','Gq'})
hold on;
s = scatter([repmat(1,1,length(Gi)),repmat(2,1,length(Gq))],[Gi Gq],'filled',...
    'jitter','on','jitterAmount',.1);
s.MarkerFaceAlpha = .5;
[p,h] = ranksum(Gi,Gq)
ylim([min([Gi Gq])-5000, max([Gi Gq])+5000])
text(1,max([Gi Gq])+1000,strcat('p =',{' '},num2str(p)),'FontSize',14)
ylabel('mean fluorescent pixel intensity');
title('63x images')
set(gca,'fontsize',15)


%% Compare fluorescent intensity Gq vs Gi (5x whole-slice images)

GiPath = 'E:\Staining Results\DREADD Gi\WEKA\For_FluoIntensity_Analysis\FluoIntensity_Images_ForMatlab';
GiFolders = dir(GiPath);
GiFolders = GiFolders(3:end);

for m = 1:length(GiFolders)
    path = strcat(GiPath,'\',GiFolders(m).name);
    images = dir(path);
    images = images(3:end);
    for i = 1:length(images)
        data = imread(strcat(path,'\',images(i).name));
        imdata_Gi(i).(GiFolders(m).name) = data; 
    end
end

allVals_Gi = []; %each C is a mouse, each R is a image mean value
for f = 1:length(fieldnames(imdata_Gi))
    mice = fieldnames(imdata_Gi);
    allIms = [];
    for i = 1:length(imdata_Gi)
        image = imdata_Gi(i).(mice{f});
        val = mean(find(image)); %finnd mean pixel value of NON-ZERO pixels
        allIms = [allIms ; val];
    end
    allVals_Gi = [allVals_Gi allIms];
end

GqPath = 'E:\Staining Results\DREAAD (Gq + LFP) Project\WEKA_Classifiers\For_FluoIntensity_Analysis\FluoIntensity_Images_ForMatlab';
GqFolders = dir(GqPath);
GqFolders = GqFolders(3:end);

for m = 1:length(GqFolders)
    path = strcat(GqPath,'\',GqFolders(m).name);
    images = dir(path);
    images = images(3:end);
    for i = 1:length(images)
        data = imread(strcat(path,'\',images(i).name));
        imdata_Gq(i).(GqFolders(m).name) = data; 
    end
end

allVals_Gq = []; %each C is a mouse, each R is a image mean value
for f = 1:length(fieldnames(imdata_Gq))
    mice = fieldnames(imdata_Gq);
    allIms = [];
    for i = 1:length(imdata_Gq)
        image = imdata_Gq(i).(mice{f});
        val = mean(find(image)); %finnd mean pixel value of NON-ZERO pixels
        allIms = [allIms ; val];
    end
    allVals_Gq = [allVals_Gq allIms];
end    

Gi = mean(allVals_Gi);
Gq = mean(allVals_Gq);

figure;
boxplot([Gi Gq],[repmat(1,1,length(Gi)),repmat(2,1,length(Gq))],'Labels',{'Gi','Gq'})
hold on;
s = scatter([repmat(1,1,length(Gi)),repmat(2,1,length(Gq))],[Gi Gq],'filled',...
    'jitter','on','jitterAmount',.1);
s.MarkerFaceAlpha = .5;
[p,h] = ranksum(Gi,Gq)
ylim([min([Gi Gq])-5e5, max([Gi Gq])+5e5])
text(1,max([Gi Gq])+1000,strcat('p =',{' '},num2str(p)),'FontSize',14)
ylabel('mean fluorescent pixel intensity');
title('5x images')
set(gca,'fontsize',15)
