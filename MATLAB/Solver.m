function [A,L] = Solver(X,Ai,Li,lambda_A,lambda_L,eta,p,n,epsilon)

% X: px(n+1)      Observed Factors
% Ai,Li: pxp         Initial Solutions
% epsilon : (typically 10^-8 to 10^-5) Stopping Criterion


% x(i) - x(i-1)
DX = zeros(p,n);
for i=1:n
    DX(:,i) = X(:,i+1)-X(:,i);
end
% Y = sum (x(i) - x(i-1))*x(i-1)^T
Y = zeros(p,p);
for i=1:n
    Y = Y + DX(:,i)*X(:,i)';
end
Y = (eta/n)*Y;
% Q = sum x(i-1)*x(i-1)^T
Q = zeros(p,p);
for i=1:n
    Q = Q + X(:,i)*X(:,i)';
end
Q = (eta^2/n)*Q;


% "p" corresponds to the "p"revious values in the iteration
A = Ai; Ap = Ai;
L = Li; Lp = Li;
t = 1;  tp = 1;

while 1
    Y_A = A + ((tp-1)/t).*(A - Ap);
    Y_L = L + ((tp-1)/t).*(L - Lp);
    
    Loss_Update = (Y_A+Y_L)*Q - Y;
    
    G_L = Y_L - Loss_Update;
    [U S V] = svd(G_L);
    for i=1:max(size(S))
        if (S(i,i)>lambda_L)
            S(i,i) = S(i,i) - lambda_L;
        else
            S(i,i)=0;
        end
    end
    Lp = L;
    L = U*S*V';
    
    G_A = Y_A - Loss_Update;
    Ap = A;
    for i=1:max(size(A))
        for j=1:max(size(A))
            if (abs(G_A(i,j))>lambda_A)
                A(i,j) = G_A(i,j) - sign(G_A(i,j))*lambda_A;
            else
                A(i,j)=0;
            end
        end
    end
    
    tp = t;
    t = (1+sqrt(4*t^2+1))/2;
    
    if ((sum(abs(A-Ap))+sum(abs(L-Lp)))/(sum(abs(A))+sum(abs(L))) < epsilon)
        break
    end
end



