% This script is to create classifiers to classify two different goup of
% patients.

clear
close all

addpath(genpath('C:\Users\xi\Box Sync\xlibrary\MATLAB'))
addpath(genpath('C:\Users\XiHang\Box Sync\xlibrary\MATLAB'))
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

%% Load negative class data
load(fileName_data_neg);
stackData_neg = stackData;
normalizedStackData_neg = normalizedStackData;
load(fileName_A_neg)
x0_neg = x0est_best;
A_neg = Aest_best;

%% Load postive class data
load(fileName_data_pos)
stackData_pos = stackData;
normalizedStackData_pos = normalizedStackData;
load(fileName_A_pos)
x0_pos = x0est_best;
A_pos = Aest_best;

%% Use the original Data for training and testing
stackData_neg = stackData_neg(:,2:end,:);
stackData_pos = stackData_pos(:,2:end,:);

mean_all = zeros(1,size(stackData_neg,2));
std_all = zeros(1,size(stackData_neg,2));
min_all = zeros(1,size(stackData_neg,2));
max_all = zeros(1,size(stackData_neg,2));

stackData_all = cat(3,stackData_neg,stackData_pos);
stackData_all_collapse = reshape(permute(stackData_all,[1,3,2]),[],size(stackData_all,2));

for k = 1:size(stackData_all,2)
    temp1 = stackData_all(:,k);
    
    mean_all(k) = mean(temp1(~isnan(temp1)));
    std_all(k) = std(temp1(~isnan(temp1)));
    min_all(k) = min(temp1(~isnan(temp1)));
    max_all(k) = max(temp1(~isnan(temp1)));
    
end

stackData = stackData_neg;
means = zeros(size(stackData,3),size(stackData,2));
stds = zeros(size(stackData,3),size(stackData,2));
mins = zeros(size(stackData,3),size(stackData,2));
maxs = zeros(size(stackData,3),size(stackData,2));
for k = 1:size(means,1)
    temp = squeeze(stackData(:,:,k));
    for m = 1:size(means,2)
        temp1 = temp(:,m);
        
        tempVal = mean(temp1(~isnan(temp1)));
        if isempty(tempVal)
            tempVal = min_all(m);
        elseif isnan(tempVal)
            tempVal = mean_all(m);
        end
        means(k,m) = tempVal;
        
        tempVal = std(temp1(~isnan(temp1)));
        if isempty(tempVal)
            tempVal = std_all(m);
        elseif isnan(tempVal)
            tempVal = std_all(m);
        end
        stds(k,m) = tempVal;
        
        tempVal = min(temp1(~isnan(temp1)));
        if isempty(tempVal)
            tempVal = min_all(m);
        elseif isnan(tempVal)
            tempVal = min_all(m);
        end
        mins(k,m) = tempVal;
        
        tempVal = max(temp1(~isnan(temp1)));
        if isempty(tempVal)
            tempVal = max_all(m);
        elseif isnan(tempVal)
            tempVal = max_all(m);
        end
        maxs(k,m) = tempVal;
    end
end

x_raw_neg = [means, stds, mins, maxs];

% positive
stackData = stackData_pos;
means = zeros(size(stackData,3),size(stackData,2));
stds = zeros(size(stackData,3),size(stackData,2));
mins = zeros(size(stackData,3),size(stackData,2));
maxs = zeros(size(stackData,3),size(stackData,2));
for k = 1:size(means,1)
    temp = squeeze(stackData(:,:,k));
    for m = 1:size(means,2)
        temp1 = temp(:,m);
        
        tempVal = mean(temp1(~isnan(temp1)));
        if isempty(tempVal)
            tempVal = min_all(m);
        elseif isnan(tempVal)
            tempVal = mean_all(m);
        end
        means(k,m) = tempVal;
        
        tempVal = std(temp1(~isnan(temp1)));
        if isempty(tempVal)
            tempVal = std_all(m);
        elseif isnan(tempVal)
            tempVal = std_all(m);
        end
        stds(k,m) = tempVal;
        
        tempVal = min(temp1(~isnan(temp1)));
        if isempty(tempVal)
            tempVal = min_all(m);
        elseif isnan(tempVal)
            tempVal = min_all(m);
        end
        mins(k,m) = tempVal;
        
        tempVal = max(temp1(~isnan(temp1)));
        if isempty(tempVal)
            tempVal = max_all(m);
        elseif isnan(tempVal)
            tempVal = max_all(m);
        end
        maxs(k,m) = tempVal;
    end
end

x_raw_pos = [means, stds, mins, maxs];

%% contruct training and testing set
label = [zeros(size(x0_neg,2),1);ones(size(x0_pos,2),1)];
A_neg_collapse = reshape(A_neg,[],size(A_neg,3));
A_pos_collapse = reshape(A_pos,[],size(A_pos,3));

perm_ind = randperm(size(A_neg,3));
x_train_neg_ind = perm_ind(1:floor(end/2));
x_test_neg_ind = perm_ind(floor(end/2)+1:end);
y_train_neg = zeros(length(x_train_neg_ind),1);
y_test_neg = zeros(length(x_test_neg_ind),1);

perm_ind = randperm(size(A_pos,3));
x_train_pos_ind = perm_ind(1:floor(end/2));
x_test_pos_ind = perm_ind(floor(end/2)+1:end);
y_train_pos = ones(length(x_train_pos_ind),1);
y_test_pos = ones(length(x_test_pos_ind),1);

x_train_neg = [x0_neg(:,x_train_neg_ind)' A_neg_collapse(:,x_train_neg_ind)'];
x_test_neg = [x0_neg(:,x_test_neg_ind)' A_neg_collapse(:,x_test_neg_ind)'];
x_train_pos = [x0_pos(:,x_train_pos_ind)' A_pos_collapse(:,x_train_pos_ind)'];
x_test_pos = [x0_pos(:,x_test_pos_ind)' A_pos_collapse(:,x_test_pos_ind)'];

% x_train_neg = [x0_neg(:,x_train_neg_ind)'];
% x_test_neg = [x0_neg(:,x_test_neg_ind)' ];
% x_train_pos = [x0_pos(:,x_train_pos_ind)' ];
% x_test_pos = [x0_pos(:,x_test_pos_ind)' ];


x_train = [x_train_neg; x_train_pos];
x_test = [x_test_neg; x_test_pos];


x_train_raw_neg = x_raw_neg(x_train_neg_ind,:);
x_test_raw_neg = x_raw_neg(x_test_neg_ind,:);
x_train_raw_pos = x_raw_pos(x_train_pos_ind,:);
x_test_raw_pos = x_raw_pos(x_test_pos_ind,:);
x_train_raw = [x_train_raw_neg; x_train_raw_pos];
x_test_raw = [x_test_raw_neg; x_test_raw_pos];

y_train = [y_train_neg; y_train_pos];
y_test = [y_test_neg; y_test_pos];

% normalize data

% x_train_normalized = x_train;
% x_test_normalized = x_test;
% x_train_raw_normalized = x_train_raw;
% x_test_raw_normalized = x_test_raw;

[x_train_normalized, setting] = ...
    normalize_lr(x_train,ones(size(x_train,1),1));
x_test_normalized = ...
    normalize_lr(x_test,ones(size(x_test,1),1),setting);
[x_train_raw_normalized, setting] = ...
    normalize_lr(x_train_raw,ones(size(x_train,1),1));
x_test_raw_normalized = ...
    normalize_lr(x_test_raw,ones(size(x_test_raw,1),1),setting);

%% Use PCA to reduce the dimensionality of the features
% xtr = x_train_normalized;
% xte = x_test_normalized;
% xtr_raw = x_train_raw_normalized;
% xte_raw = x_test_raw_normalized;
% 
% out_tr = pcadm(xtr);
% x_train_normalized = out_tr.rdcddata;
% out_te = pcadm(xte,[],out_tr);
% x_test_normalized = out_te.rdcddata;
% 
% out_tr_raw = pcadm(xtr_raw);
% x_train_raw_normalized = out_tr_raw.rdcddata;
% out_te_raw = pcadm(xte_raw,[],out_tr_raw);
% x_test_raw_normalized = out_te_raw.rdcddata;


%% Learn classifier

theta = logiRegr(x_train_normalized,y_train);
y_pred = logiEval(theta,x_test_normalized);
perf= binary_classification_performance_evaluation(y_test,y_pred);

theta_raw = logiRegr(x_train_raw_normalized,y_train);
y_pred_raw = logiEval(theta_raw,x_test_raw_normalized);
perf_raw= binary_classification_performance_evaluation(y_test,y_pred_raw);

svmModel = fitcsvm(x_train_normalized,y_train,'Standardize',true,'KernelFunction','rbf',...
    'KernelScale','auto');
[~,y_pred_svm] = predict(svmModel,x_test_normalized);
perf_svm = binary_classification_performance_evaluation(y_test,y_pred_svm(:,2));

svmModel = fitcsvm(x_train_raw_normalized,y_train,'Standardize',true,'KernelFunction','rbf',...
    'KernelScale','auto');
[~, y_pred_svm] = predict(svmModel,x_test_raw_normalized);
perf_svm_raw= binary_classification_performance_evaluation(y_test,y_pred_svm(:,2));




