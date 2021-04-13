function [INT, IDX, f_FFT, FFT, t_FFT, v_FFT, swTimes] = sleepDetection_zscore_multFiles(data, sf, names, figOutput)
%% 
%Inputs: 
% data: cell array, each entry is from CSV file: matrix with each column corresponding to a recording
% channel, C1 = time points, C2 = LFP, C3 = EMG, C4 = optoswitch (see lines
% 12-16)
% sf: Sampling frequency
% name: cell array containing each name of dataset, used to label figures
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

numTrials = length(data);
idx = [];
allLFP = [];
for t = 1:numTrials
    currLFP = data{t}(:,2);
    allLFP = [allLFP ; currLFP];
    idx = [idx t*ones(1, length(currLFP))];
end
lfpZ = zscore(allLFP);
    

%% Calculate spectogram
params.Fs = sf;
params.tapers = [5 9]; %[1.5 1];  % 
params.fpass = [.1 100]; %frequencies to analyze
params.err = [2 0.05];
params.trialave = 0;
params.pad = 0;
movingwin =  [10 5]; %[3 .1]; %[10 10];  %change?? [win winstep]

FFT = {}; t_FFT = {}; f_FFT = {};
for t = 1:numTrials
    currLFP = lfpZ(idx == t);
    [currFFT,currt,currf,~]=mtspecgramc(currLFP,movingwin,params);
    FFT{t} = currFFT;
    t_FFT{t} = currt;
    f_FFT{t} = currf;
end

%% Calculate locomotion
notMoving = {}; v_FFT = {};
for t = 1:numTrials
    opto = round(data{t}(:,4));
    v_FFT_curr = []; %locomotion: 1 if moving, 0 if not

    window = movingwin(1);   %in sec
    winstep = movingwin(2);  %in sec: Each bin will be spaced apart by winstep
    for i = 1:length(t_FFT{t})
        m=t_FFT{t}(i) - (winstep/2); %t_FFT are the midpoints of the bins, but we want to index from the beginning of a bin
        mm = m*(sf); %original LFP sampling freq. 
        Optowindow = mm + (winstep*sf-1);
        diffWin = diff(opto(round(mm:Optowindow))); 
        if ~isempty(find(diffWin, 1)) 
            v_FFT_curr(i) = 1; %moving
        else 
            v_FFT_curr(i) = 0; % not moving 
        end
    end

    currNotMoving = (v_FFT_curr' == 0); %1 if animal is not moving
    notMoving{t} = currNotMoving;
    v_FFT{t} = v_FFT_curr;
end

%% Calculate EMG

window = movingwin(1);   %in sec
winstep = movingwin(1)-movingwin(2);  %in sec

%using all trials, take the zscore of all EMG values 
allEMG = []; idx = [];
for t = 1:numTrials
    currEMG = data{t}(:,3);
    allEMG = [allEMG ; currEMG];
    idx = [idx t*ones(1, length(currEMG))];
end
emgZ = abs(zscore(allEMG)); 

emgAll = {};
for t = 1:numTrials
    emgAll{t} = emgZ(idx == t);
end

t_EMG = {}; emg = {};
for t = 1:numTrials
    t_FFT_curr = t_FFT{t};
%     currEMG = abs(zscore(data{t}(:,3)));
    currEMG = emgAll{t};
    %EMG must be same length as t_FFT for later - the following takes the mean for each t_FFT time bin:
    meanEMG = [];
    for i = 1:length(t_FFT_curr)
        m = t_FFT_curr(i) - (winstep/2); %t_FFT are the midpoints of the bins, but we want to index from the beginning of a bin
        mm = m*(sf); %original LFP sampling freq. 
        EMGwindow = mm + (winstep*sf);
        meanEMG = [meanEMG mean(currEMG(round(mm:EMGwindow)))];
    end
    % emg = abs(meanEMG'); %take the absolute value of this so that we can pick a threshold more easily
    emg{t} = meanEMG;
end

emgThresh = 0.4; %means EMG was + or - 0.25 SD from mean of EMG (aka baseline)

%% Isolate NREM sleep with delta power zscore

swsFreq = [0.5 4]; %each row is a freq range to be used to pickout SWS
otherFreq = [8 20]; %compar slow wave to this other freq range for ratio

%find corresponding idx in f
swRatioAll = []; idx = [];
for t = 1:numTrials
    [~,f_first] = min(abs(f_FFT{t}-swsFreq(1)));
    [~,f_last] = min(abs(f_FFT{t}-swsFreq(2)));
    sw_f = [f_first:f_last];

    [~,f_first] = min(abs(f_FFT{t}-otherFreq(1)));
    [~,f_last] = min(abs(f_FFT{t}-otherFreq(2)));
    other_f = [f_first:f_last];

    swPower = mean(log(FFT{t}(:,sw_f)),2);
    otherPower = mean(log(FFT{t}(:,other_f)),2);
    swRatio = otherPower ./ swPower;
    swRatioAll = [swRatioAll ; swRatio];
    idx = [idx t*ones(1, length(swRatio))];
end 
swRatioZ = zscore(swRatioAll);

swRatio = {};
for t = 1:numTrials
    swRatio{t} = smooth(swRatioZ(idx == t),15);
end

swthresh = 0; %threshold of zscored slow wave power to find sleep periods
swthresh_win_sec = [10 5]; %[win step] in sec, window to examine mean zscore
swthresh_win = swthresh_win_sec/movingwin(2);

NREMemgThresh = 5;

NREMtimes = {}; swTimes = {};
for t = 1:numTrials    
    %For moving window:
    NREMtimesCurr = []; %binary for NREM or Other
    swTimesCurr = []; %corresponds to NREMtimes, has the mid-bin time point
    i=1;
    while i < length(swRatio{t}) - (swthresh_win(1)-1)
        st = i;
        fin = i + swthresh_win(1)-1;
        swTimesCurr = [swTimesCurr (t_FFT{t}(st)+t_FFT{t}(fin))/2]; %new mid-bin time point

        if mean(swRatio{t}(st:fin))>swthresh && isempty(find(notMoving{t}(st:fin)-1)) %&& mean(emg{t}(st:fin))<NREMemgThresh %criteria for NREM sleep! 
            NREMtimesCurr = [NREMtimesCurr 1];
        else
            NREMtimesCurr = [NREMtimesCurr 0];
        end

        i = i + swthresh_win(2);
    end
    NREMtimes{t} = NREMtimesCurr; 
    swTimes{t} = swTimesCurr;
end

%% Find theta power for REM sleep

thFreq = [6 10]; %each row is a freq range to be used to pickout SWS

thPowerAll = []; idx = [];
for t = 1:numTrials    
    %find corresponding idx in f
    first = thFreq(1);
    last = thFreq(2);
    [~,f_first] = min(abs(f_FFT{t}-first));
    [~,f_last] = min(abs(f_FFT{t}-last));
    th_f = f_first:f_last;
    
    thPower = mean(log(FFT{t}(:,th_f)),2); 
    thPowerAll = [thPowerAll ; thPower];
    idx = [idx t*ones(1, length(thPower))];
end
thRatioZ = zscore(thPowerAll);

thPower = {};
for t = 1:numTrials
    thPower{t} = smooth(thRatioZ(idx == t),15);
end

thThresh = 0.25; %threshold of zscored slow wave power to find sleep periods
thThresh_win_sec = [10 5]; %[win step] in sec, window to examine mean zscore MUST BE SAME AS NREM IF YOU WANT TO USE BUSZAKI CODE
thThresh_win = thThresh_win_sec/movingwin(2);

REMtimes = {}; thTimes = {};
for t = 1:numTrials
    %For moving window:
    REMtimesCurr = [];
    thTimesCurr = [];
    i=1;
    while i < (length(thPower{t}) - (thThresh_win(1)-1))
        st = i;
        fin = i + thThresh_win(1)-1;
        thTimesCurr = [thTimesCurr (t_FFT{t}(st)+t_FFT{t}(fin))/2]; %new mid-bin time point

        if mean(thPower{t}(st:fin))>thThresh && isempty(find(notMoving{t}(st:fin)-1)) && mean(emg{t}(st:fin))<emgThresh && NREMtimes{t}(i) == 0
            REMtimesCurr = [REMtimesCurr 1];
        else
            REMtimesCurr = [REMtimesCurr 0];
        end

        i = i + thThresh_win(2);
    end
    REMtimes{t} = REMtimesCurr;
    thTimes{t} = thTimesCurr;
end

%% Plot

halfWin = movingwin(2)/2;
halfWin_NREM = swthresh_win_sec(2)/2;
halfWin_REM = thThresh_win_sec(2)/2;

if figOutput
    for t = 1:numTrials
        times = data{t}(:,1);
        lfp = data{t}(:,2);
        
        figure;
        hold on;
        ha(1) = subplot(3,1,1);
            hold on;
            title(names{t})
            imagesc([t_FFT{t}(1) t_FFT{t}(end)], [f_FFT{t}(1) f_FFT{t}(end)], log(FFT{t}(:,10:end)'))
            ylim([0 20])
            caxis([-10 0]);
        ha(2) = subplot(3,1,2);
            plot(times/1000,lfp)
        ha(3)= subplot(3,1,3);
            hold on;
            plot(t_FFT{t}, thPower{t})
            plot(t_FFT{t}, swRatio{t})
            plot(t_FFT{t}, v_FFT{t})
            legend({'theta', 'slow wave', 'locomotion'})
            for n=1:length(swTimes{t})
                if NREMtimes{t}(n)==1
                    s = swTimes{t}(n) - halfWin_NREM;
                    fin = swTimes{t}(n) + halfWin_NREM;
                    h = fill([s s fin fin], [-2 4 4 -2], [0 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
                end
            end
            for n = 1:length(thTimes{t})
                if REMtimes{t}(n) == 1
                    s = thTimes{t}(n) - halfWin_REM;
                    fin = thTimes{t}(n) + halfWin_REM;
                    h = fill([s s fin fin], [-2 4 4 -2], [0 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
                end
            end
            linkaxes(ha,'x')
    end
end

%% Produce outputs similar to buszaki code
% ONLY WORKS IF REM AND NREM MOVING WINDOW ARE THE SAME AND SO SWTIMES AND
% REMTIMES ARE THE SAME

IDX = {}; INT = {};
for t = 1:numTrials
    IDX{t} = NREMtimes{t} + 2*REMtimes{t} +1;
    INT{t} = IDXtoINT(IDX{t},3); %Currently has start and end points in swTimes space 

    for i = 1:length(INT{t})
        onsets = INT{t}{i}(:,1); %in swTimes space
        onsets = swTimes{t}(onsets)'; %now in sec but corresponds to mid-point of bin
        onsets = onsets - halfWin_NREM; %adjusted to be real onset now

        offsets = INT{t}{i}(:,2);
        offsets = swTimes{t}(offsets)';
        offsets = offsets + halfWin_NREM; %adjusted offsets too

        INT{t}{i} = [onsets offsets];
    end

    wake = find(IDX{t} == 1);
    sws = find(IDX{t} == 2);
    rem = find(IDX{t} == 3);

end


