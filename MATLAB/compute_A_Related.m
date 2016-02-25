function [A2n, IphA2n, dIphA2ndA] = compute_A_Related (A,h,N)

A2n = cell(N+1,1); % because of the case of A^0
IphA2n = cell(N+1,1); % because of the case of IphA2n^0
dIphA2ndA = cell(N+1,1);


A2n{1} = eye(size(A,1));
for k = 1:N
    A2n{k+1} = A^k;
end


IphA2n{1} = eye(size(A,1));
for k = 1:N
    IphA2n{k+1} = eye(size(A,1));
    for j = 1:k
%         IphA2n{k+1} = IphA2n{k+1} + nchoosek(k,j)*h^j*A2n{(j)+1};
        IphA2n{k+1} = IphA2n{k+1} + nchoosek_times_h2k(k,j,h)*A2n{(j)+1};
    end
end

dAndA = cell(N+1,1);
dAndA{1} = zeros(size(A,1)^2);

for k = 1:N
    dAndA{k+1} = dAn_dA(A,k);
end

dIphA2ndA{1} = zeros(size(A,1)^2);
for ni = 1:N
    dIphA2ndA{ni+1} = zeros(size(A,1)^2);
    for k = 1:ni
%         dIphA2ndA{ni+1} = dIphA2ndA{ni+1} + nchoosek(ni,k)*h^k*dAndA{k+1};
        dIphA2ndA{ni+1} = dIphA2ndA{ni+1} + nchoosek_times_h2k(ni,k,h)*dAndA{k+1};
    end
end

function out = nchoosek_times_h2k(k,j,h)
% replace nchoosek(k,j)*h^j to nchoosek_times_h2k(k,j,h)

a = k:-1:(k-j+1);
b = 1:j;
c = a./b *h;
out = prod(c); 





