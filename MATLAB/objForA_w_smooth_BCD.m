function [f,g] = objForA_w_smooth_BCD(vecA,Aprev,n,x0,xt,h,lambda)

nd = sqrt(length(vecA));

A = reshape(vecA,nd,nd);
I = eye(nd);

measuredx0 = xt(:,1);

xt0 = x0;
xt = xt(:,2:end);

m = length(n);
N = n(end);
kronx0I = kron(xt0',I);

% B0 is a vector to record whether this are any missing in the initial
% measurement
B0 = double(diag(~isnan(measuredx0)));
measuredx0(isnan(measuredx0)) = 0;

% B is a cell array which record whether there are missings in each time
% point of the measurements, each B matrix in a cell is a diagonal matrix
% with zero's and ones; one indicate the measurement is valid, and 0 is
% invalid
B = cell(m,1);
for i = 1:m
    B{i} = double(diag(~isnan(xt(:,i))));
    xt(isnan(xt(:,i)),i) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%
% Calculate utilities %
%%%%%%%%%%%%%%%%%%%%%%%

[~, IphA2n, dIphA2ndA] = compute_A_Related (A,h,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function value %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Error Term -- 
termError = 0;
for i = 1:m
    termError = termError ...
        + (IphA2n{n(i)+1}*xt0 - xt(:,i))'*B{i}*(IphA2n{n(i)+1}*xt0 - xt(:,i));
end
termError = 1/2*termError;

% -- regularization term --
termRegularization = lambda.reg/2*(vecA - reshape(Aprev,[],1))'*(vecA - reshape(Aprev,[],1));

% -- smoothing term --
termSmooth = 0;
for k = 1:N
    termSmooth = termSmooth + ...
        (IphA2n{k+1}*xt0-IphA2n{k}*xt0)'*(IphA2n{k+1}*xt0-IphA2n{k}*xt0);
end
termSmooth = termSmooth*lambda.smooth/2;

% -- initial term error --
termInitial = 0.5*(x0-measuredx0)'*B0*(x0-measuredx0)*lambda.initial;

% -- Objective function value --
f = termError + termRegularization + termSmooth + termInitial;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Differential of Objective function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Error Term -- 
dtermError = 0;

for i = 1:m
    dtermError = dtermError + ...
        (IphA2n{n(i)+1}*xt0-xt(:,i))'*B{i}*kronx0I*dIphA2ndA{n(i)+1};
end
dtermError = dtermError';

% -- regularization term --
dtermRegularization = lambda.reg*(vecA - reshape(Aprev,[],1));

% -- smoothing term --
dtermSmooth = 0;
for k = 1:N
    dtermSmooth = dtermSmooth + ...
        (xt0'*IphA2n{k+1}'*kronx0I - xt0'*IphA2n{k}'*kronx0I)*...
        (dIphA2ndA{k+1}-dIphA2ndA{k});
%         (IphA2n{k}*xt0)'*kronx0I*dIphA2ndA{k+1} + ...
%         (IphA2n{k+1}*xt0)'*kronx0I*dIphA2ndA{k};
end
dtermSmooth = lambda.smooth*dtermSmooth';

% -- First Differential of Objective function
g = dtermError + dtermRegularization + dtermSmooth;

% disp(['f = ' num2str(f)]);