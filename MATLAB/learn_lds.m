clear
close all

load('../data/gse54514_top10_processed.mat')

nSub = size(normalizedStackData,3);
nBiomarker = size(normalizedStackData,2)-1;
subID = 1:nSub;

chosen_subID = subID(chosen_ind);

Aest_best = zeros(nBiomarker,nBiomarker,length(subID));
x0est_best = zeros(nBiomarker,length(subID));
dt = 0.2;
param.lambda_reg = 1e-4;
param.lambda_smooth = 0;
param.lambda_initial = 1;
param.dt = dt;

n_chosen = length(chosen_subID);
cA = cell(n_chosen,1);
cfval = cell(n_chosen,1);
cx0 = cell(n_chosen,1);

for k = chosen_subID
    tStart = tic;
    fprintf(['Case #' num2str(k)])
%     disp(['worker#' num2str(k) ' is working'])
    tData = squeeze(normalizedStackData(:,:,k));
    tData = tData';
    nanTp = isnan(tData(2,:));  % use 2 because if one variable is missing other variables will be missing too
    tData(:,nanTp) = [];
%     ts = tData(1,:);

    Aest = zeros(size(tData,1)-1);
    N = 3;
    A = zeros(nBiomarker,nBiomarker,N);
    fval = zeros(N,1);
    x0 = zeros(nBiomarker,N);
    inputTData = tData;
    for m = 1:N
        t1 = tic;
        try
            [A(:,:,m), x0(:,m),fval(m)] = optDepMtx_BCD(inputTData,Aest,param);
        catch
            A(:,:,m) = inf(nBiomarker, nBiomarker);
            x0(:,m) = inf(nBiomarker,1);
            fval(m) = inf;
        end
        t2 = toc(t1);
        disp(['Window # ' num2str(m) ' takes ' num2str(t2) ' fVal = ' num2str(fval(m))])
    end
    cA{k} = A;
    cx0{k} = x0;
    cfval{k} = fval;
    lapse = toc(tStart);
    fprintf([' used time ' num2str(lapse) '\n']);
end