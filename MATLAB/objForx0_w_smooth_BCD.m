function [f,g] = objForx0_w_smooth_BCD(xt0,Aprev,IphA2n,xt,n,h,lambda)

nd = length(xt0);

measuredx0 = xt(:,1);
xt = xt(:,2:end);
m = length(n);
N = n(end);

vecA = reshape((IphA2n{2}-eye(nd))/h,[],1); % IphA2n{1} = eye(nd);

B0 = double(diag(~isnan(measuredx0)));
measuredx0(isnan(measuredx0)) = 0;
B = cell(m,1);
for i = 1:m
    B{i} = double(diag(~isnan(xt(:,i))));
    xt(isnan(xt(:,i)),i) = 0;
end

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
termInitial = 0.5*(xt0-measuredx0)'*B0*(xt0-measuredx0)*lambda.initial;

% -- Objective function value --
f = termError + termRegularization + termSmooth + termInitial;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Differential of Objective function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Error Term -- 
dtermError = 0;

for i = 1:m
    dtermError = dtermError + ...
        (IphA2n{n(i)+1}*xt0-xt(:,i))'*B{i}*IphA2n{n(i)+1};
end
dtermError = dtermError';

% -- smoothing term --
dtermSmooth = 0;
for k = 1:N
    dtermSmooth = dtermSmooth + ...
        xt0'*(IphA2n{k+1}-IphA2n{k})'*(IphA2n{k+1}-IphA2n{k}); 
end
dtermSmooth = lambda.smooth*dtermSmooth';

% -- term initial error
dtermInitial = B0*(xt0-measuredx0)*lambda.initial;

% -- First Differential of Objective function
g = dtermError + dtermSmooth + dtermInitial;






