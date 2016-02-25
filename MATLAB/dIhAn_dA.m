function out = dIhAn_dA (A,h,n)

p = size(A,1);

out = zeros(p^2);

for k = 1:n
    
    out = out + nchoosek(n,k) * h^k * dAn_dA(A,k);
    
end