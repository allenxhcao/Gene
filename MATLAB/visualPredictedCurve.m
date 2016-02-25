% This script is used to visualize the prediction curve after learning the
% interaction network

clear
close all

groupName_neg = 'infection';
groupName_pos = 'shock_pos';

selection_crit = '_age60after_1weeks_2subsamples';

fileName_data_neg = ['data/labevents_sepsis_' ...
    groupName_neg ...
    '_hadm_id_adult_random_500_Stack_cbc' selection_crit];
fileName_A_neg = ['result/A_individualslabevents_sepsis_'...
    groupName_neg ...
    '_hadm_id_adult_random_500_Stack_cbc' selection_crit '_3'];

fileName_data_pos = ['data/labevents_sepsis_' ...
    groupName_pos ...
    '_hadm_id_adult_random_500_Stack_cbc' selection_crit];
fileName_A_pos = ['result/A_individualslabevents_sepsis_'...
    groupName_pos ...
    '_hadm_id_adult_random_500_Stack_cbc' selection_crit '_3'];

load(fileName_data_neg);
load(fileName_A_neg);

Aest_best_all = Aest_best;
x0est_best_all = x0est_best;

normalizedStackData(:,redundantBiomarkerID+1,:)=[];
biomarkerName(redundantBiomarkerID+1) = [];
allTs = squeeze(normalizedStackData(:,1,:));
allTs  = allTs(:);
tsLength = max(allTs(~isnan(allTs)));
nT = tsLength/dt;

pred_all_neg = zeros(nT+1,length(biomarkerName),length(subID));

for iSub = subID
    
    Aest_best = squeeze(Aest_best_all(:,:,iSub));
    x0est_best = squeeze(x0est_best_all(:,iSub));
    tData_neg = squeeze(normalizedStackData(:,:,iSub));
    ts = tData_neg(1,:);
    
    tData_neg = tData_neg';
    nanTp = isnan(tData_neg(1,:));
    tData_neg(:,nanTp) = [];
    
    pred = dlds(Aest_best,x0est_best,nT,dt);
    pred(2:end,:) = max(min(pred(2:end,:),1),0);
    
    pred_all_neg(:,:,iSub) = pred';
    
    visual_true_pred(tData_neg,pred);
    
end
%%
load(fileName_data_pos);
load(fileName_A_pos);

Aest_best_all = Aest_best;
x0est_best_all = x0est_best;

normalizedStackData(:,redundantBiomarkerID+1,:)=[];
biomarkerName(redundantBiomarkerID+1) = [];
nMarker = size(Aest_best_all,1);
allTs = squeeze(normalizedStackData(:,1,:));
allTs  = allTs(:);
tsLength = max(allTs(~isnan(allTs)));
nT = tsLength/dt;

pred_all_pos = zeros(nT+1,length(biomarkerName),length(subID));

for iSub = subID
    
    Aest_best = squeeze(Aest_best_all(:,:,iSub));
    x0est_best = squeeze(x0est_best_all(:,iSub));
    tData_pos = squeeze(normalizedStackData(:,:,iSub));
    ts = tData_pos(1,:);
    
    tData_pos = tData_pos';
    nanTp = isnan(tData_pos(1,:));
    tData_pos(:,nanTp) = [];
    
    pred = dlds(Aest_best,x0est_best,nT,dt);
    pred(2:end,:) = max(min(pred(2:end,:),1),0);
    
    pred_all_pos(:,:,iSub) = pred';
    
end
 
pred_all_neg_mean = mean(pred_all_neg,3)';
pred_all_pos_mean = mean(pred_all_pos,3)';

% pred_all_neg_mean = max(pred_all_neg,[],3)';
% pred_all_pos_mean = max(pred_all_pos,[],3)';

% pred_all_neg_mean = min(pred_all_neg,[],3)';
% pred_all_pos_mean = min(pred_all_pos,[],3)';

%% Plot the mean of all
figure('Units','pixels','Position',[0 0 1680 1080]);
for k = 2:length(biomarkerName)
    subplot(ceil(nMarker/5),5,k-1)
    hold on
    for iSub = 1:size(pred_all_neg,3)
        pred = squeeze(pred_all_neg(:,:,iSub))';
        plot(pred(1,:),pred(k,:),'Color',[0,1,0,0.1],'LineWidth',2);
    end
    
    for iSub = 1:size(pred_all_pos,3)
        pred = squeeze(pred_all_pos(:,:,iSub))';
        plot(pred(1,:),pred(k,:),'Color',[1,0,0,0.1],'LineWidth',2);
    end
    
    plot(pred_all_neg_mean(1,:),pred_all_neg_mean(k,:),'Color',[0,1,0],'LineWidth',5)
    plot(pred_all_pos_mean(1,:),pred_all_pos_mean(k,:),'Color',[1,0,0],'LineWidth',5)
    title([biomarkerName{k}]);
    xlabel('Hour');
    ylim([0 1])
    xlim([0 tsLength])
    grid on
    hold off
end

