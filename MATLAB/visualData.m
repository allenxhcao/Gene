clear
close all

groupName_neg = 'infection';
groupName_pos = 'shock_pos';

selection_crit = '_age60after_1weeks_2subsamples';

fileName_data_neg = ['data/labevents_sepsis_' ...
    groupName_neg ...
    '_hadm_id_adult_random_500_Stack_cbc' selection_crit];
fileName_data_pos = ['data/labevents_sepsis_' ...
    groupName_pos ...
    '_hadm_id_adult_random_500_Stack_cbc' selection_crit];

load(fileName_data_neg);

nMarker = size(stackData,2)-1;
dt = 1;
stackData_neg = stackData;
allTs = squeeze(normalizedStackData(:,1,:));
allTs  = allTs(:);
tsLength = max(allTs(~isnan(allTs)));
nT = tsLength/dt;

%%
load(fileName_data_pos);
stackData_pos = stackData;


%%
figure('Units','pixels','Position',[0 0 1680 1080]);
for k = 2:length(biomarkerName)
    subplot(ceil(nMarker/5),5,k-1)
    hold on
    for iSub = 1:size(stackData_neg,3)
        data = squeeze(stackData_neg(:,:,iSub))';
        plot(data(1,:),data(k,:),'o','Color',[0,1,0,0.1],'LineWidth',2);
    end
    
    for iSub = 1:size(stackData_pos,3)
        data = squeeze(stackData_pos(:,:,iSub))';
        plot(data(1,:),data(k,:),'o','Color',[1,0,0,0.1],'LineWidth',2);
    end
    
%     plot(pred_all_neg_mean(1,:),pred_all_neg_mean(k,:),'Color',[0,1,0],'LineWidth',5)
%     plot(pred_all_pos_mean(1,:),pred_all_pos_mean(k,:),'Color',[1,0,0],'LineWidth',5)
    title([biomarkerName{k}]);
    xlabel('Hour');
    xlim([0 tsLength])
    grid on
    hold off
end
    
%     Apos = Aest_best.*(Aest_best>=0);
%     Aneg = abs(Aest_best.*(Aest_best<=0));
%     [t1,~] = graphminspantree(sparse(Apos));
%     t1 = full(t1);
%     [t2,~] = graphminspantree(sparse(Aneg));
%     t2 = full(t2);
%     A = t1 + -1*t2 + Aest_best.*eye(size(Aest_best));
%     [bg,~] = drawGraph(A,biomarkerName(2:end));
%     view(bg)
