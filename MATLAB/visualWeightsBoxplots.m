% This script is to visualize the box plot of the weights and initial
% conditions of the results from learning time series data

clear
% close all


addpath(genpath('C:\Users\xi\Box Sync\xlibrary\MATLAB'))
addpath(genpath('C:\Users\XiHang\Box Sync\xlibrary\MATLAB'))
groupName_neg = 'infection';
groupName_pos = 'shock_pos';

selection_crit = '_age60after_1weeks_2subsamples_2';

fileName_A_neg = ['result/A_individualslabevents_sepsis_'...
    groupName_neg ...
    '_hadm_id_adult_random_500_Stack_cbc' selection_crit];

fileName_A_pos = ['result/A_individualslabevents_sepsis_'...
    groupName_pos ...
    '_hadm_id_adult_random_500_Stack_cbc' selection_crit];

%% Load negative class data
load(fileName_A_neg)
x0_neg = x0est_best;
A_neg = Aest_best;

%% Load postive class data
load(fileName_A_pos)
x0_pos = x0est_best;
A_pos = Aest_best;

nMarker = size(x0_neg,1);
%% box plot of initial condition
% figure('Units','pixels','Position',[100 0 500 200*nMarker])
% for k = 1:nMarker
%     subplot(nMarker,1,k)
%     data_neg = x0_neg(k,:)';
%     data_pos = x0_pos(k,:)' ;
%     data = [data_neg; data_pos];
%     group = [zeros(size(data_neg)); ones(size(data_pos))];
%     boxplot(data,group,'label',{'Negative','Positive'});
%     
% end

%% box plot of network weights
figure('Units','pixels','Position',[100 0 200*nMarker 200*nMarker])
count = 0;
for k = 1:nMarker
    for m = 1:nMarker
        count = count + 1;
        subplot(nMarker,nMarker,count)
        data_neg = squeeze(A_neg(k,m,:));
        data_pos = squeeze(A_pos(k,m,:));
        data = [data_neg; data_pos];
        group = [zeros(size(data_neg)); ones(size(data_pos))];
        boxplot(data,group,'label',{'Negative','Positive'});
        ylim([-1 1])
    end
end





























