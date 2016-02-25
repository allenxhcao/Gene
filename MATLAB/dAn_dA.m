function out = dAn_dA (A,n)

%% approach one

out = zeros(size(A,1)^2);

for j = 1:n
    out = out + kron((A')^(n-j),A^(j-1));
end
