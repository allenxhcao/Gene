function visual_true_pred(tData, pred)

nMarker = size(tData,1)-1;

figure('Units','pixels','Position',[200 0 200 nMarker*200]);
for k = 1:nMarker
    subplot(nMarker,1,k)
    hold on
    plot(tData(1,:),tData(k+1,:),'ro');
    plot(pred(1,:),pred(k+1,:),'-b');
    hold off
end
close