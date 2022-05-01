%子空间迭代法求特征值
function [w_m,w] = SSI(q,K,M,N_w)
err_tol = 1e-8;
n = size(K,1);
X1 = zeros(n,q);
X1(:,1)=ones(n,1);
for i = 2:q
    X1(i-1,i) = 1;
end
X1p = X1;
Y1p = M*X1p;
X2 = K\Y1p;
Y2 = M*X2;
K_bar = X2'*Y1p;
M_bar = X2'*Y2;
[phi_bar2,gama_bar2] = eig(M_bar\K_bar);
gama_bar1 = eye(q);
k=1;

while max(abs(diag(gama_bar2)-diag(gama_bar1))./diag(gama_bar2))>err_tol
    k=k+1;
    if k>1e4 
        fprintf('超出迭代次数限制\n');
    break;
    end
    Y2p = Y2*phi_bar2;
    Y1p = Y2p;
    X2 = K\Y1p;
    Y2 = M*X2;
    K_bar = X2'*Y1p;
    M_bar = X2'*Y2;
    
    gama_bar1 = gama_bar2;
    [phi_bar2,gama_bar2] = eig(M_bar\K_bar);
    for i=1:q  
        phi_bar2(:,i)=phi_bar2(:,i)/sqrt(phi_bar2(:,i)'*M_bar*phi_bar2(:,i));
    end
end
lamda = gama_bar2;
[lamda,~] = sort(diag(lamda));
w = (lamda(N_w))^0.5;
w_m = (max(lamda))^0.5;


end
