function [out, lfp_delta, lfp_gamma, isflipped] = detectSOvDelta_121319(lfp,Fs,sleepTimes,figOn)
%% [out] = detect_slow_oscillations(lfp,Fs)
%   identifies times of down states during sleep block
%       1. averages activity across channels
%       2. filters average lfp signal for delta activity (.1-4Hz)
%       3. finds positive-to-negative zero crossings (zc) during sleep
%       4. for each pos-to-neg zc
%           4.1. get the prev peak (down state)
%           4.2. get the prev negative-to-positive zc
%           4.3. get the next negative-to-positive zc
%       5. returns slow oscillations that fit the following criteria
%           5.1. time btw neg-to-pos zcs is btw .8 and 2s
%           5.2. prev peak exceeds a threshold of 90% of all peaks
%       
%   Inputs:
%       lfp - lfp is a matrix (samples x channels) of lfp data from sleep
%       Fs - sample rate of lfp 
% 
%   Optional Name/Value Pairs:
%       'sleepTimes' - a logical vector 
%                   - (1-sleep,0-awake) in seconds
%       'artifact_idx' - a logical vector 
%                   - (1-artifact,0-fine) corresponding to that session's lfp
%       'bad_channels' - vector of channels that are ignored (default = [])
%       'figOn' - [0|1]
%       'sleep_classify' - [0|1], only for classified sleep (sleep_idx)
%       'mnl_parm' - [peak-thr trough-thr dur-min dur-max], manual parameter setting
% 
%   Outputs:
%       out
%           .down_states - time of down states
%           .zc - time of pos-to-neg down state
%           .peak - peak of z-scored delta band lfp signal
%           .trough - trough of z-scored delta band lfp signal

%% INPUTS
MinPeakHeight = 90;
mnl_parm=[85 40 .15 .5];

%% ORGANIZE LFP
isflipped = 0;
% sizing
dt = 1/Fs;
[N,num_ch] = size(lfp);
time = ((1:N)/Fs)' - dt;

lfp = double(zscore(lfp));

%% CREATE SLEEPTIMES THAT MATCHES LFP LENGTH

sleepIdx = zeros(length(lfp),1);

for k = find(sleepTimes) %for each sec of sleep
    sleepIdx((k*Fs+1):(k*Fs+Fs)) = 1;
end

%% FILTER LFP FOR DELTA
if length(mnl_parm)==4, fpass = [.1,4];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif length(mnl_parm)==6, fpass=mnl_parm(5:6); end
lfp_delta = filter_delta(lfp,Fs,fpass);

%% POS-TO-NEG ZERO CROSSINGS
idx = round(Fs):(N-round(Fs));
ptnzc = round(Fs)-1 + find(lfp_delta(idx)>=0 & lfp_delta(idx+1)<0 & sleepIdx(idx)==1);

%% GET STATS ON ALL POS-TO-NEG ZERO CROSSINGS
prev_peak_idx = zeros(size(ptnzc));
prev_peak     = zeros(size(ptnzc));
next_trough   = zeros(size(ptnzc));
dur           = zeros(size(ptnzc));
isasleep      = zeros(size(ptnzc));
prev_ntpzc = zeros(size(ptnzc));
next_ntpzc = zeros(size(ptnzc));
for i=2:length(ptnzc)-1
    prev_ntpzc(i)          = ptnzc(i-1)-1 + find(lfp_delta(ptnzc(i-1):ptnzc(i  ))<0,1,'last');
    next_ntpzc(i)          = ptnzc(i)  -1 + find(lfp_delta(ptnzc(i  ):ptnzc(i+1))<0,1,'last');
    [prev_peak(i),idx]  = max(lfp_delta(ptnzc(i-1):ptnzc(i)));
    prev_peak_idx(i)    = ptnzc(i-1)-1 + idx;
    [next_trough(i),idx]= min(lfp_delta(ptnzc(i):next_ntpzc(i)));
    next_trough_idx(i)  = ptnzc(i)-1 + idx;
%     dur(i)              = time(next_ntpzc - prev_ntpzc);
    dur(i)              = time(next_trough_idx(i) - prev_peak_idx(i));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    isasleep(i)         = sleepIdx(prev_peak_idx(i)) & sleepIdx(next_trough_idx(i));
end

%% USE GAMMA TO DETERMINE UP VS DOWN STATES
thresh_up = prctile(prev_peak,mnl_parm(1));%85%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresh_dwn = prctile(next_trough,mnl_parm(2));%40%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresh = [thresh_up,thresh_dwn];

idx = isasleep & prev_peak>=thresh_up & next_trough<thresh_dwn & dur>mnl_parm(3) & dur<mnl_parm(4);%SO
idx2 = isasleep & prev_peak<thresh_up & next_trough<thresh_dwn & dur>mnl_parm(3)-.05; %delta %03/09/2019

idxTest = idx + idx2;

% Bandpass filter 
fpass = [80 100];

n = round((3/fpass(1)) * Fs); %set order equal to 3 cycles of low cutoff frequency

b = designfilt('bandpassfir', 'FilterOrder', n, 'CutoffFrequency1', fpass(1), 'CutoffFrequency2', fpass(2), 'SampleRate', Fs);

lfp_gamma = filtfilt(b,double(lfp));
lfp_gamma = abs(zscore(lfp_gamma)); %removed abs 06/19/20 (why was it there?!)
%%%%%%%%%%%%%%%%%%%%

peak_g = []; trough_g = [];
for i = [find(idxTest)]'
    peak_g = [peak_g mean(lfp_gamma(prev_ntpzc(i):ptnzc(i)))];
    trough_g = [trough_g mean(lfp_gamma(ptnzc(i):next_ntpzc(i)))];
end

if mean(peak_g) > mean(trough_g) %if gamma during peak (down state) is more than trough (up state), then lfp has to be flipped
    isflipped = 1; %mark as 1 to report that lfp was flipped
    lfp = -1*(lfp);
    lfp_delta = -1*(lfp_delta);
    lfp_gamma = -1*lfp_gamma;
    %POS-TO-NEG ZERO CROSSINGS
    idx = round(Fs):(N-round(Fs));
    ptnzc = round(Fs)-1 + find(lfp_delta(idx)>=0 & lfp_delta(idx+1)<0 & sleepIdx(idx)==1);
    %GET STATS ON ALL POS-TO-NEG ZERO CROSSINGS
    prev_peak_idx = zeros(size(ptnzc));
    prev_peak     = zeros(size(ptnzc));
    next_trough   = zeros(size(ptnzc));
    dur           = zeros(size(ptnzc));
    isasleep      = zeros(size(ptnzc));
    prev_ntpzc = zeros(size(ptnzc));
    next_ntpzc = zeros(size(ptnzc));
    for i=2:length(ptnzc)-1
        prev_ntpzc(i)          = ptnzc(i-1)-1 + find(lfp_delta(ptnzc(i-1):ptnzc(i  ))<0,1,'last');
        next_ntpzc(i)          = ptnzc(i)  -1 + find(lfp_delta(ptnzc(i  ):ptnzc(i+1))<0,1,'last');
        [prev_peak(i),idx]  = max(lfp_delta(ptnzc(i-1):ptnzc(i)));
        prev_peak_idx(i)    = ptnzc(i-1)-1 + idx;
        [next_trough(i),idx]= min(lfp_delta(ptnzc(i):next_ntpzc(i)));
        next_trough_idx(i)  = ptnzc(i)-1 + idx;
        dur(i)              = time(next_trough_idx(i) - prev_peak_idx(i));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        isasleep(i)         = sleepIdx(prev_peak_idx(i)) & sleepIdx(next_trough_idx(i));
    end
    temp = peak_g; %flip these
    peak_g = trough_g;
    trough_g = temp;
end   

thresh_up = prctile(prev_peak,mnl_parm(1));%85%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresh_dwn = prctile(next_trough,mnl_parm(2));%40%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresh = [thresh_up,thresh_dwn];

%% APPLY CRITERIA FOR SLOW OSCILLATIONS

% idx = isasleep & prev_peak>=thresh_up & next_trough<thresh_dwn & dur>.3 & dur<1;
idx = isasleep & prev_peak>=thresh_up & next_trough<thresh_dwn & dur>mnl_parm(3) & dur<mnl_parm(4);%.15%.5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%delta %03/09/2019
idx2 = isasleep & prev_peak<thresh_up & next_trough<thresh_dwn & dur>mnl_parm(3)-.05;

%% organize outputs
%SO
out.down_states = time(prev_peak_idx(idx));
out.up_states = time(next_trough_idx(idx));
out.zc = time(ptnzc(idx));
out.peaks = prev_peak(idx);
out.troughs = next_trough(idx);
%delta %03/09/2019
out.down_states_d = time(prev_peak_idx(idx2));
out.up_states_d = time(next_trough_idx(idx2));
out.peaks_d = prev_peak(idx2);
out.troughs_d = next_trough(idx2);
out.zc_d = time(ptnzc(idx2));

%% plot
if figOn,
    figure;
    cc = get(gca,'ColorOrder');
        
    clf
    set(gcf,'Name',sprintf('Slow_Oscillations'))
    set(gcf,'Position',[680,620,980,360])
    
    subplot(2,3,1), hold on
    histogram(peak_g);
    histogram(trough_g);
    legend({'DOWN state gamma','UP state gamma'})
    
    t = 0:.1:10;
    subplot(2,3,4), hold on
    isi = diff(out.down_states);
    histogram(isi,0:.2:10,'Normalization','pdf')
    lam = lognfit(isi);
    plot(t,lognpdf(t,lam(1),lam(2)),'LineWidth',2)
    [mx,idx] = max(lognpdf(t,lam(1),lam(2)));
    plot(t(idx),mx,'k*','MarkerSize',10)
    text(t(idx)+1,mx,sprintf('Peak @ %.2fHz',1/t(idx)))
    xlabel('sec')
    ylabel('Autocorrelation')
    
    ax(1)=subplot(2,3,2:3); hold on
    plot(time,lfp,'Color',cc(1,:))
    plot(time(find(sleepIdx)),repmat(1*max(lfp),sum(sleepIdx),1),'.','Color',cc(2,:))
    
    ax(2)=subplot(2,3,5:6); hold on
    plot(time,lfp_delta,'Color',cc(1,:))
    plot(out.down_states,out.peaks,'k*')
    plot(out.down_states_d, .15*ones(1,length(out.down_states_d)),'r*')
    hline(thresh,'k-')
    hline(0,'k--')
    
    linkaxes(ax,'x')
end

end % detect down states

function lfp_delta = filter_delta(lfp,Fs,fpass)
    % filter for delta
%     fpass = [.1,4];

    % highpass
    [b,a] = butter(2,fpass(1)/(Fs/2),'high');
    lfp_delta = filtfilt(b,a,lfp);

    % lowpass
    [b,a] = butter(4,fpass(2)/(Fs/2),'low');%original from Daniel
%     [b,a] = butter(2,fpass(2)/(Fs/2),'low');%JK 04/18/19
    lfp_delta = filtfilt(b,a,lfp_delta);
end

