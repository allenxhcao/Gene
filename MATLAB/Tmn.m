function tmn = Tmn(m,n)

vec = (1:m*n)';

mtx = reshape (vec,m,[]);

mtxTran = mtx';

t = eye(m*n);

tmn = t(reshape(mtxTran,[],1),:);

