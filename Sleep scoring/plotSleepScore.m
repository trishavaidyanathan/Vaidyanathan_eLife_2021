function  [sws,wake,rem] = plotSleepScore(IDX,swTimes,iszscore,iscolorbar)

%% input zscore = 1 if IDX was obtained using sleepscore_zscore code
%% input zscore = 0 if IDX was obtained using buszaki code
%% input colobar = 1 if you want colorbar, 0 if not

wake = find(IDX == 1);
sws = find(IDX == 2);
rem = find(IDX == 3);
    
if iszscore
%With my zscore code:
    scatter(swTimes(wake), -IDX(wake), 'b', 'filled')
    hold on;
    scatter(swTimes(sws), -IDX(sws), 'r', 'filled')
    hold on;
    scatter(swTimes(rem), -IDX(rem),'g', 'filled')
    ylim([-4 0])
    xlim([0 3600])
    set(gca,'YTick',[-3:-1])
    set(gca,'YTickLabel',{'REM','SWS','Wake/MA'})
%     set(gca,'XTickLabel',{})
    if iscolorbar 
        colorbar;
    end
else
   %With buzsaki code:
    scatter(swTimes(wake), -IDX(wake), 'b', 'filled')
    hold on;
    scatter(swTimes(sws), -IDX(sws), 'r', 'filled')
    hold on;
    scatter(swTimes(rem), -IDX(rem),'g', 'filled')
    ylim([-4 0])
    xlim([0 3000])
    set(gca,'YTick',[-3:-1])
    set(gca,'YTickLabel',{'REM','SWS','Wake/MA'})
    %set(gca,'XTickLabel',{})
    if iscolorbar
        iscolorbar;
    end
end
