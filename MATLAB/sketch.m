clear
close all

fName = 'labevents_sepsis_shock_pos_hadm_id_adult_random_500_Stack_cbc';
load(['data/' fName]);
redundantBiomarkerID = [];
normalizedStackData(:,redundantBiomarkerID+1,:)=[];
biomarkerName(redundantBiomarkerID+1) = [];

nSub = size(normalizedStackData,3);
nBiomaker = size(normalizedStackData,2)-1;
subID = 1:nSub;

Aest_best = zeros(nBiomaker,nBiomaker,length(subID));
x0est_best = zeros(nBiomaker,length(subID));

param.lambda = 1e-4;
param.dt = 1;

for k = 1:length(subID)
    tStart = tic;
    fprintf(['Case #' subName{k} '\n'])
%     disp(['worker#' num2str(k) ' is working'])
    tData = normalizedStackData(:,:,subID(k));
    tData = tData';
    dt = 1;
    nanTp = isnan(tData(1,:));
    tData(:,nanTp) = [];
    ts = tData(1,:);
    % incr = 10; % increment of time units (dt) for each window
    % wSize = ceil(max(ts)/param.dt); % window size in terms of time units (dt)
    % 
    % firstTs = tData(1,1); % the time of the first time  point
    % lastTs = tData(1,end); % the time of last time point
    % 
    % tsStart = firstTs:incr*dt:lastTs-wSize*dt;
    % tsEnd = tsStart + wSize*dt - dt;

    Aest = zeros(size(tData,1)-1);
    N = 3;
    cA = cell(N,1);
    cfval = zeros(N,1);
    cx0 = cell(N,1);
    inputTData = tData;
    n = cumsum(round(diff(inputTData(1,:))/1)); 
    N = n(end);
end
