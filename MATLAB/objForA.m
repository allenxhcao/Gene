function [f,g] = objForA(vecA,Aprev,n,xt,h,lambda)

nd = sqrt(length(vecA));

A = reshape(vecA,nd,nd);
I = eye(nd);
xt0 = xt(:,1);
xt = xt(:,2:end);

m = length(n);

term1 = 0;

for i = 1:m
    
    term1 = term1 ...
        + ((I + h*A)^n(i)*xt0 - xt(:,i))'*((I + h*A)^n(i)*xt0 - xt(:,i));
    
end

term1 = 1/2*term1;

term2 = lambda/2*(vecA - reshape(Aprev,[],1))'*(vecA - reshape(Aprev,[],1));

f = term1 + term2;

dJ1 = 0;

for i = 1:m
    
    dJ1 = dJ1 + ((I + h*A)^n(i)*xt0 - xt(:,i))'*kron(xt0',I)*dIhAn_dA(A,h,n(i));
    
end

dJ1 = dJ1';

dJ2 = lambda*(vecA - reshape(Aprev,[],1));

g = dJ1 + dJ2;

% disp(['f = ' num2str(f)]);