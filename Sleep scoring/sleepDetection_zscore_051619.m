function [INT, IDX, f, FFT, t_FFT, v_FFT, swTimes] = sleepDetection_zscore_051619(data, sf, name, figOutput)
%% 
%Inputs: 
% data: From CSV file, matrix with each column corresponding to a recording
% channel so C1 = time points, C2 = LFP, C3 = EMG, C4 = optoswitch (see lines
% 12-16)
% sf: Sampling frequency
% name: name of dataset, used to label figures
% figOutput: 1 if you want figures plotted, 0 if not

%Outputs: 
% INT: Cell array containing onsets and offsets of wake (INT{1}), nrem
% (INT{2}), rem (INT{3}) in seconds
% IDX: Each entry corresponds to swTimes, contains a 1 for wake, 2 for
% nrem, 3 for rem
% f: list of frequencies in spectogram output FFT
% FFT: spectogram output
% t_FFT: list of time mid-points in sec that correspond to spectogram
% output FFT
% v_FFT: binarized locomotion (1 is moving, 0 is stationary), corresponds
% to time bins in t_FFT
% swTimes: mid-point of time bins in sec that corresponds to IDX


%% Load data

times = data(:,1); %in ms
lfp = data(:,2);
emg = data(:,3);
opto = round(data(:,4));

%% Calculate spectogram
params.Fs = sf;
params.tapers = [5 9];  
params.fpass = [.1 100]; %frequencies to analyze
params.err = [2 0.05];
params.trialave = 0;
params.pad = 0;
movingwin =  [10 5]; 

[FFT,t_FFT,f,~]=mtspecgramc(data(:,2),movingwin,params);

%% Calculate locomotion
v_FFT = []; %locomotion: 1 if moving, 0 if not

window = movingwin(1);   %in sec
winstep = movingwin(2);  %in sec: Each bin will be spaced apart by winstep
for i = 1:length(t_FFT)
    m=t_FFT(i) - (winstep/2); %t_FFT are the midpoints of the bins, but we want to index from the beginning of a bin
    mm = m*(sf); %original LFP sampling freq. 
    Optowindow = mm + (winstep*sf-1);
    diffWin = diff(opto(mm:Optowindow)); 
    if ~isempty(find(diffWin, 1)) 
        v_FFT(i) = 1; %moving
    else 
        v_FFT(i) = 0; % not moving 
    end
end

notMoving = (v_FFT' == 0); %1 if animal is not moving

%% Calculate EMG
window = movingwin(1);   %in sec
winstep = movingwin(1)-movingwin(2);  %in sec
meanEMG = [];
EMGz = (abs(zscore(data(:,3)))); %EMG amplitude
%EMG must be same length as t_FFT for later - the following takes the mean for each t_FFT time bin:
for i = 1:length(t_FFT)
    m=t_FFT(i) - (winstep/2); %t_FFT are the midpoints of the bins, but we want to index from the beginning of a bin
    mm = m*(sf); %original LFP sampling freq. 
    EMGwindow = mm + (winstep*sf);
    meanEMG = [meanEMG mean(EMGz(mm:EMGwindow))];
end
emg = meanEMG;
t_EMG = t_FFT;

emgThresh = 0; %means EMG was + or - 0.25 SD from mean of EMG (aka baseline)

%% Isolate NREM sleep with delta power zscore

swsFreq = [0.5 4]; %each row is a freq range to be used to pickout SWS
otherFreq = [8 20]; %compare slow wave to this other freq range for ratio

%find corresponding idx in f
sw_f = [];
for i = 1:size(swsFreq,1)
    first = swsFreq(i,1);
    last = swsFreq(i,2);
    [~,f_first] = min(abs(f-first));
    [~,f_last] = min(abs(f-last));
    sw_f = [sw_f f_first:f_last];
end

[~,f_first] = min(abs(f-otherFreq(1)));
[~,f_last] = min(abs(f-otherFreq(2)));
other_f = [f_first:f_last];

swPower = mean(log(FFT(:,sw_f)),2);
otherPower = mean(log(FFT(:,other_f)),2);
swRatio = otherPower ./ swPower;
swRatio = smooth(zscore((swRatio)),15);

swthresh = 0.5; %threshold of zscored slow wave power to find sleep periods
swthresh_win_sec = [10 5]; %[win step] in sec, window to examine mean zscore
swthresh_win = swthresh_win_sec/movingwin(2);

%For moving window:
NREMtimes = []; %binary for NREM or Other
swTimes = []; %corresponds to NREMtimes, has the mid-bin time point
i=1;
while i < length(swRatio) - (swthresh_win(1)-1)
    st = i;
    fin = i + swthresh_win(1)-1;
    swTimes = [swTimes (t_FFT(st)+t_FFT(fin))/2]; %new mid-bin time point
    
    if mean(swRatio(st:fin))>swthresh && isempty(find(notMoving(st:fin)-1))
        NREMtimes = [NREMtimes 1];
    else
        NREMtimes = [NREMtimes 0];
    end
    
    i = i + swthresh_win(2);
end

%% Find theta power for REM sleep

thFreq = [6 10]; %each row is a freq range to be used to pickout SWS

%find corresponding idx in f
th_f = [];
for i = 1:size(thFreq,1)
    first = thFreq(i,1);
    last = thFreq(i,2);
    [~,f_first] = min(abs(f-first));
    [~,f_last] = min(abs(f-last));
    th_f = [th_f f_first:f_last];
end

thPower = smooth(zscore(mean(log(FFT(:,th_f)),2)),15); %picking out .5-3 hz (~delta) and 9-12 hz (spindles)

thThresh = 0.25; %threshold of zscored slow wave power to find sleep periods
thThresh_win_sec = [10 5]; %[win step] in sec, window to examine mean zscore MUST BE SAME AS NREM IF YOU WANT TO USE BUSZAKI CODE
thThresh_win = thThresh_win_sec/movingwin(2);

%For moving window:
REMtimes = [];
thTimes = [];
i=1;
while i < length(thPower) - (thThresh_win(1)-1)
    st = i;
    fin = i + thThresh_win(1)-1;
    thTimes = [thTimes (t_FFT(st)+t_FFT(fin))/2]; %new mid-bin time point
    
    if mean(thPower(st:fin))>thThresh && isempty(find(notMoving(st:fin)-1)) && mean(emg(st:fin))<emgThresh && NREMtimes(i) == 0
        REMtimes = [REMtimes 1];
    else
        REMtimes = [REMtimes 0];
    end
    
    i = i + thThresh_win(2);
end


%% Plot

    halfWin = movingwin(2)/2;
    halfWin_NREM = swthresh_win_sec(2)/2;
    halfWin_REM = thThresh_win_sec(2)/2;
if figOutput
    figure;
    ha(1) = subplot(3,1,1);
    title(name)
    hold on;
    imagesc([t_FFT(1) t_FFT(end)], [f(1) f(end)], log(FFT(:,10:end)'))
    ylim([0 20])
    ha(2) = subplot(3,1,2);
    plot(times/1000,lfp)
    ha(3)= subplot(3,1,3);
    plot(t_FFT,thPower)
    hold on;
    plot(t_FFT,swRatio)
    hold on;
    plot(t_FFT, v_FFT)
    legend({'theta', 'slow wave', 'locomotion'})
    hold on
    for n=1:length(swTimes)
        if NREMtimes(n)==1
            s = swTimes(n) - halfWin_NREM;
            fin = swTimes(n) + halfWin_NREM;
            h = fill([s s fin fin], [-2 4 4 -2], [0 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            hold on;
        end
    end
    hold on;
    for n = 1:length(thTimes)
        if REMtimes(n) == 1
            s = thTimes(n) - halfWin_REM;
            fin = thTimes(n) + halfWin_REM;
            h = fill([s s fin fin], [-2 4 4 -2], [0 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
            hold on;
        end
    end
    linkaxes(ha,'x')
end

%% Produce outputs similar to buszaki code
% ONLY WORKS IF REM AND NREM MOVING WINDOW ARE THE SAME AND SO SWTIMES AND
% REMTIMES ARE THE SAME

IDX = NREMtimes + 2*REMtimes +1;
INT = IDXtoINT(IDX,3); %Currently has start and end points in swTimes space 

for i = 1:length(INT)
    onsets = INT{i}(:,1); %in swTimes space
    onsets = swTimes(onsets)'; % in sec but corresponds to mid-point of bin
    onsets = onsets - halfWin_NREM; %adjusted to be real onset 
    
    offsets = INT{i}(:,2);
    offsets = swTimes(offsets)';
    offsets = offsets + halfWin_NREM; %adjusted offsets too
    
    INT{i} = [onsets offsets];
end

wake = find(IDX == 1);
sws = find(IDX == 2);
rem = find(IDX == 3);
