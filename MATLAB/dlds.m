function out = dlds(A,x0,t,dt)

p = length(x0); % number of dimensions

ts = 0:dt:t;
nt = length(ts); % number of time points to generate

out = zeros(p,nt);
out(:,1) = x0;

for k = 2:nt
    
    temp = (eye(p) + dt*A)*out(:,k-1);
%     temp(temp<0) = 0;
    out(:,k) = temp;
        
end

out = [ts;out];



