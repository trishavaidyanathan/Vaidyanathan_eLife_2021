function [pre, post, Spect, meanTraces, T, F] = CaETSpect(datastruct, resstruct, events, window, islfp, name)

%events = cell array, same size of datastruct, each element corresponds to a row in datastruct and 
%contains evt indeces to be evaulated

%Example of how I used it:
% % datastruct = allDataCyto;
% % resstruct = allResCyto;
% % name = 'cyto sleep'
% % events = {};
% % for i = 1:length(datastruct)
% %     events{i} = [datastruct(i).sleepEvtIdx{:}]
% % end
% % [Sepctmed_norm, T, F] = ETSpect(allDataCyto, allResCyto, events, name)
% % 


Spect = [];

movingwin = [2,0.1]; %[2,0.1]; %moving window for spectrogram

%LFP params structure
Params.Fs = 1000;
Params.fpass = [0.1 100]; 
Params.pad =-1;
Params.tapers = [2.5 2];
Params.trialave = 0;

freqs = [0.5,4;10,20;20,50;70,100];
f = 1;

pre = []; post = []; Spect = [];
for w=1:length(datastruct) %take a structure of files with astrocyte and LFP data  
    if islfp
        lfp = datastruct(w).lfp;
    else
        lfp = datastruct(w).eeg;
    end
    if ~isempty(lfp) %if there is lfp data    
        event_times = resstruct(w).fts.curve.tBegin(events{w})*datastruct(w).frameperiod; %get an array of onsets and convert to sec
        event_times = event_times(event_times > window(1) & event_times < (length(lfp)/Params.Fs)-window(2)); %make sure window doesn't exceed matrix

        [S,T,F] = mtspecgramtrigc(lfp,event_times,window,movingwin,Params); %event triggered spectrogram
        S = single(S); 
        if ~isempty(S)
            Spect =cat(3,Spect,S); %concatenate
        end

        clear lfp S event_times
        disp([num2str(w) '/' num2str(length(datastruct))]) %track prog
    end
end

Spectmed_norm = Spect ./ median(Spect); %normalzie in time
Fidx = F>=freqs(f,1) & F<=freqs(f,2);
Spect = Spectmed_norm(:,Fidx,:); %select out desired frequencies
meanTraces = squeeze(mean(Spect,2)); %average across freqs, each C is the ETA trace for an event

Spectmed_norm = median(Spectmed_norm,3);

figure; plot_matrix(Spectmed_norm,T,F,'nl'); %plot spectrogram of relative power
caxis([.92 1.06])
line(window,[0 Params.fpass(2)],'linewidth',2,'color','k') %plot spectrogram of relative power
title(name);

half = (length(T)-1)/2;
[~,twosec] = min(abs(T-2));
pre = [mean(meanTraces(half-twosec:half,:))]; %two sec before
post = [mean(meanTraces(half+1:half+1+twosec,:))];

end
