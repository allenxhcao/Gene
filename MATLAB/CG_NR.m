function [optx,niter] = CG_NR(x,f,fp,fpp)

if nargin < 3
    fp = @(x)(gradest(f,x)');
    fpp = @(x)(hessian(f,x)');
elseif nargin < 4
    fpp = @(x)(hessian(f,x)');
end


imax = 20e4;
jmax = 10e4;
n = length(x);
t1 = 1e-8;
t2 = 1e-4;

i = 0;
k = 0;
r = -fp(x);
d = r;
deltaNew = r'*r;
deltaOld = inf;
% fOld = inf;
% fNew = f(x);

while (i < imax) && (abs(deltaOld - deltaNew) > abs(t1*deltaNew) && deltaNew > t1)
    j = 0;
    deltaD = d'*d;
    alpha = -(fp(x)'*d)/(d'*fpp(x)*d);
    x = x + alpha * d;
    j = j + 1;
    while j < jmax && alpha^2*deltaD > t2^2
        alpha = -(fp(x)'*d)/(d'*fpp(x)*d);
        x = x + alpha * d;
        j = j + 1;
    end
    r = -fp(x);
    deltaOld = deltaNew;
    deltaNew = r'*r;
    beta = deltaNew/deltaOld;
    d = r + beta*d;
    k = k + 1;
    if k == n || r'*d<=0
        d = r;
        k = 0;
    end
    i = i+1;
%     fOld = fNew;
%     fNew = f(x);
end
optx = x;
niter = i;
    
        

    


