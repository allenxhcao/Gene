% tData: multi-variate time series data. Each row is a measurement
% Apre: dependend matrix of previous window
% param: parameters

function [A,fVal] = optDepMtx(tData, Aprev, param)

lambda = param.lambda;
h = param.dt;
% 
% c1 = param(1);
% c2 = param(2);


[nRow,~] = size(tData); % nt: number of measurements; nd: number of dimensions

% nt = nCol - 1;  % number of time points except the initial one
nd = nRow - 1;  % number of variables (dimensions)

% find all time units between all consecutive examples 
n = cumsum(round(diff(tData(1,:))/h)); 
xt = tData(2:end,:);    % time data

objFun = @(vecA)objForA_w_smooth(vecA,Aprev,n,xt,h,lambda);

options = optimoptions('fminunc',...
        'GradObj','on',...
        'Algorithm','quasi-newton',...
        'Display','Iter',...
        'MaxIter',2000,...
        'TolFun',1e-5);
% options = optimoptions('fminunc',...
%         'GradObj','on',...
%         'Algorithm','quasi-newton',...
%         'Display','iter-detailed',...
%         'HessUpdate','dfp',...
%         'MaxIter',4000);
    

if sum(Aprev(:))~=0
    vecA0 = reshape(Aprev,[],1);
else
    if exist('referenceANoTExist.mat','file')
%     if exist('referenceA.mat','file')
        load('referenceA');
        vecA0 = reshape(referenceA,[],1);
    else
        vecA0 = (rand(nd^2,1)-0.5-reshape(eye(nd^2,1),[],1))*0.1;
        
    end
end

[vecA,fVal] = fminunc(objFun,vecA0,options);

A = reshape(vecA,nd,[]);








% % find A(d,:) of the matrix
% for d = 1:nd
%     
%     xt = tData(2:end,d);
%     x = tData(1:end-1,:);
%     
%     f = @(a) 1/2*(x*a-xt)'*(x*a-xt) + c1/2*((a - Aprev(:,d))'*(a - Aprev(:,d)))...
%         + c2/2*(a'*a);
%     fp = @(a) x'*x*a - x'*xt + c1 * (a - Aprev(:,d)) + c2 * a;
%     fpp = @(a) x'*x + c1 * eye(nd) + c2 * eye(nd);
%     ff = {f;fp;fpp};
%     options = optimoptions('fminunc',...
%         'GradObj','on',...
%         'Algorithm','trust-region',...
%         'Hessian','on');
%     
%     if(sum(Aprev(:,d))~=0)
%         a0 = Aprev(:,d);
%     else
%         a0 = rand(nd,1);
%     end
%     [opta,~] = CG_NR(a0,f,fp,fpp);
% %     opta = fminunc(ff,Apre(:,d)+rand(nd,1),options);
%     A(:,d) = opta;
%     
% end
    