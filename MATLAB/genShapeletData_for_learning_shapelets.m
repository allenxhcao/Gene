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
    '_hadm_id_adult_random_500_Stack_cbc' selection_crit];

fileName_data_pos = ['data/labevents_sepsis_' ...
    groupName_pos ...
    '_hadm_id_adult_random_500_Stack_cbc' selection_crit];
fileName_A_pos = ['result/A_individualslabevents_sepsis_'...
    groupName_pos ...
    '_hadm_id_adult_random_500_Stack_cbc' selection_crit];

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
n_neg = length(subID);
for iSub = subID
    
    Aest_best = squeeze(Aest_best_all(:,:,iSub));
    x0est_best = squeeze(x0est_best_all(:,iSub));
    tData = squeeze(normalizedStackData(:,:,iSub));
    ts = tData(1,:);
    
    tData = tData';
    nanTp = isnan(tData(1,:));
    tData(:,nanTp) = [];
    
    pred = dlds(Aest_best,x0est_best,nT,dt);
    pred(2:end,:) = max(min(pred(2:end,:),10),-1);
    
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
n_pos = length(subID);
for iSub = subID
    
    Aest_best = squeeze(Aest_best_all(:,:,iSub));
    x0est_best = squeeze(x0est_best_all(:,iSub));
    tData = squeeze(normalizedStackData(:,:,iSub));
    ts = tData(1,:);
    
    tData = tData';
    nanTp = isnan(tData(1,:));
    tData(:,nanTp) = [];
    
    pred = dlds(Aest_best,x0est_best,nT,dt);
    pred(2:end,:) = max(min(pred(2:end,:),10),-1);
    
    pred_all_pos(:,:,iSub) = pred';
    
end

%% Generate matrix with desired format, and sub-sample
pred_all_neg(2:2:end,:,:) = [];
pred_all_pos(2:2:end,:,:) = [];

pred_all_neg(:,1,:) = [];
pred_all_pos(:,1,:) = [];

% convert the multivariate time series to single variate
% a row is a time series

way_to_convert = 'difference_concatenate';
switch way_to_convert
    case 'concatenate'
        pred_all_neg = reshape(pred_all_neg,[],n_neg)';
        pred_all_pos = reshape(pred_all_pos,[],n_pos)';

    case 'mean'
        pred_all_neg = squeeze(mean(pred_all_neg,2))';
        pred_all_pos = squeeze(mean(pred_all_pos,2))';
    
    case 'difference_concatenate'
        pred_all_neg = diff(reshape(pred_all_neg,[],n_neg))';
        pred_all_pos = diff(reshape(pred_all_pos,[],n_pos))';
end

pred_all_neg = [1*ones(n_neg,1) pred_all_neg];
pred_all_pos = [2*ones(n_pos,1) pred_all_pos];

rng(1)
ind_neg = randperm(n_neg);
ind_pos = randperm(n_pos);

ind_train_neg = ind_neg(1:round(n_neg/2));
ind_test_neg = ind_neg(round(n_neg/2)+1:end);
ind_train_pos = ind_pos(1:round(n_pos/2));
ind_test_pos = ind_pos(round(n_pos/2)+1:end);

pred_train = [pred_all_neg(ind_train_neg,:); pred_all_pos(ind_train_pos,:)];
pred_test = [pred_all_neg(ind_test_neg,:); pred_all_pos(ind_test_pos,:)];

pred_train = pred_train(randperm(size(pred_train,1)),:);
pred_test = pred_test(randperm(size(pred_test,1)),:);
%% Write text files
folderName = ['MIMICIII_' groupName_neg '_AS_NEG_VS_' groupName_pos '_AS_POS_' selection_crit '_' way_to_convert];
if ~exist(folderName,'dir')
    mkdir(folderName);
end

textFileName_train = [folderName '/' folderName '_TRAIN'];
textFileName_test = [folderName '/' folderName '_TEST'];

dlmwrite(textFileName_train,pred_train,'delimiter',',');
dlmwrite(textFileName_test,pred_test,'delimiter',',');

