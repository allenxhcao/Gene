clear
close all

groupName_neg = 'normal';
groupName_pos = 'shock_pos';

fileName_data_neg = ['data/labevents_sepsis_' ...
    groupName_neg ...
    '_hadm_id_adult_random_500_Stack_cbc'];
fileName_A_neg = ['result/A_individualslabevents_sepsis_'...
    groupName_neg ...
    '_hadm_id_adult_random_500_Stack_cbc.mat'];

fileName_data_pos = ['data/labevents_sepsis_' ...
    groupName_pos ...
    '_hadm_id_adult_random_500_Stack_cbc'];
fileName_A_pos = ['result/A_individualslabevents_sepsis_'...
    groupName_pos ...
    '_hadm_id_adult_random_500_Stack_cbc.mat'];

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
    tData = squeeze(normalizedStackData(:,:,iSub));
    ts = tData(1,:);
    
    tData = tData';
    nanTp = isnan(tData(1,:));
    tData(:,nanTp) = [];
    
    pred = dlds(Aest_best,x0est_best,nT,dt);
    pred(2:end,:) = max(min(pred(2:end,:),1),0);
    
    pred_all_neg(:,:,iSub) = pred';
    
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
    tData = squeeze(normalizedStackData(:,:,iSub));
    ts = tData(1,:);
    
    tData = tData';
    nanTp = isnan(tData(1,:));
    tData(:,nanTp) = [];
    
    pred = dlds(Aest_best,x0est_best,nT,dt);
    pred(2:end,:) = max(min(pred(2:end,:),1),0);
    
    pred_all_pos(:,:,iSub) = pred';
    
end

%% Generate matrix with desired format, and sub-sample
pred_all_neg(2:2:end,:,:) = [];
pred_all_pos(2:2:end,:,:) = [];

pred_all_neg(:,1,:) = [];
size_neg = size(pred_all_neg);
pred_all_neg = cat(1,-1*ones([1 size_neg([2,3])]),pred_all_neg);

pred_all_pos(:,1,:) = [];
size_pos = size(pred_all_pos);
pred_all_pos = cat(1,ones([1 size_pos([2,3])]),pred_all_pos);

pred = cat(3,pred_all_pos, pred_all_neg);
pred = permute(pred,[3,2,1]);
pred = reshape(pred,[],size(pred,3));

%% Write text files
textFileName = ['data/' groupName_neg '_AS_NEG_VS_' groupName_pos '_AS_POS.txt'];
dlmwrite(textFileName,pred,'delimiter',' ');

