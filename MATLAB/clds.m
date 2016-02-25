function out = clds (A,x0,t,dt)

if nargin < 4
    dt = t/100;
end

ts = dt:dt:t;

expAdt = expm(A*dt); 
out = [x0 zeros(length(x0),length(ts))];

for k = 2:size(out,2)
    temp = expAdt*out(:,k-1);
    out(:,k) = double(temp>0).*temp;
end

out = [0 ts;out];
