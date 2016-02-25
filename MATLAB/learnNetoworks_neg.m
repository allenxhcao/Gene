clear
close all

fName = 'labevents_sepsis_shock_neg_hadm_id_adult_random_500_Stack_cbc';
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

parfor k = 1:length(subID)
    tStart = tic;
    fprintf(['Case #' subName{k}])
%     disp(['worker#' num2str(k) ' is working'])
    tData = normalizedStackData(:,:,subID(k));
    tData = tData';
    dt = 1;
    nanTp = isnan(tData(1,:));
    tData(:,nanTp) = [];
    ts = tData(1,:);

    Aest = zeros(size(tData,1)-1);
    N = 3;
    cA = cell(N,1);
    cfval = zeros(N,1);
    cx0 = cell(N,1);
    inputTData = tData;
    for m = 1:N
        t1 = tic;
        try
            [cA{m}, cx0{m},cfval(m)] = optDepMtx_BCD(inputTData,Aest,param);
        catch
            cA{m} = inf(nBiomarker, nBiomarker);
            cx0{m} = inf(nBiomarker,1);
            cfval(m) = inf;
        end
        t2 = toc(t1);
%         disp(['Window # ' num2str(m) ' takes ' num2str(t2) ' fVal = ' num2str(cfval(m))])
    end
    [~,I] = min(cfval);
    Aest_best(:,:,k) = cA{I(1)};
    x0est_best(:,k) = cx0{I(1)};
    lapse = toc(tStart);
    fprintf([' used time ' num2str(lapse) '\n']);
end

save(['result/A_individuals' fName],...
    'Aest_best',...
    'x0est_best',...
    'dt',...
    'subID',...
    'redundantBiomarkerID')