%% %%%%% FIGURE 3 FIGURE SUPPLEMENT 3 %%%%%%

%% Data sets used in this figure
% PCA Gi Data
pathPCA = 'E:\Analysis\MAT files\MAT Files for Paper 2020\DATA FOR DRYAD\PCA\Gi DREADD'; %change this to your path on your computer

load(strcat(pathPCA,'\PCA_Output_Gi'));

%% Fig3FigSupp3 A 

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
'spatial overlap size corrected'}; 

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

%% Remaining figures in Fig3FigSupp3 were done in Python