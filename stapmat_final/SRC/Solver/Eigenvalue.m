%* *****************************************************************
%* - Function of STAPMAT in Eigenvalue phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To solve frequecny and eigenvalue                           *
%*                                                                 *
%* - Call procedures:                                              *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     Zixiong Lin, 2022.05.23                                     *
%*                                                                 *
%* *****************************************************************

function Eigenvalue()

global cdata;
global sdata;

NEQ = sdata.NEQ;

NUM_EIG = 10; % 求解前NUM_GIE阶的特征值
if NUM_EIG > NEQ
    NUM_EIG = NEQ;
end

sdata.FREQUENCY = zeros(NUM_EIG, 1, 'double');
sdata.PHI = zeros(NEQ, NUM_EIG, 'double');

SPSTIFF = Matrix2Sparse(sdata.STIFF);
SPMASS = Matrix2Sparse(sdata.MASS);

[V,D] = eigs(full(SPSTIFF),eye(36));
[V,D] = eig(full(SPSTIFF),full(SPMASS));
tol = 1e-6; % 收敛判据
% [x,lam,k] = inverse1(full(SPSTIFF),full(SPMASS),tol);
% [lambda,phi] = inverse(SPSTIFF,SPMASS,NUM_EIG,tol); % 特征值lambda和特征向量phi
[x,lam,k] = inverse1(full(SPSTIFF),full(SPMASS),tol);
[lambda,phi] = inverse(full(SPSTIFF),full(SPMASS),NUM_EIG,tol); % 特征值lambda和特征向量phi

end

% ----------------------- Functions -----------------------------------

% Convert the stiff vector to a sparse stiff matrix
function SPMatrix = Matrix2Sparse(A)

global sdata;
% A = sdata.STIFF;
MAXA = sdata.MAXA;
NEQ = sdata.NEQ;
NWK = sdata.NWK;
IIndex = zeros(NWK*2-NEQ, 1);
JIndex = IIndex;
STIFF = IIndex;

NUM = 1;
NUMC = 0;
for N = 1:NEQ
    KU = MAXA(N + 1) - MAXA(N);
    for L = 1:KU
        IIndex(NUM) = N;
        JIndex(NUM) = N - L + 1;
        STIFF(NUM) = A(NUM);
        NUM = NUM + 1;
        if (L == 1)
            NUMC = NUMC + 1;
            continue;
        end
        SYMN = NUM-1 - NUMC + NWK;
        IIndex(SYMN) = N - L + 1;
        JIndex(SYMN) = N;
        STIFF(SYMN) = A(NUM-1);
    end
end

SPMatrix = sparse(IIndex, JIndex, STIFF, NEQ, NEQ);
end


% 求特征值和特征向量
% 输入刚度阵K，质量阵M，求解阶数num，收敛残差tol
% 输出前num阶特征值lambda，特征向量phi
function [lambda,phi] = inverse(K,M,num,tol)

n = size(K,1);
lambda = zeros(num,1);
phi = zeros(n,num);

if abs(det(K)) < 1e-6
    alpha = 1; % 如果K不满秩则移轴
else
    alpha = 0;
end
Kbar = K + alpha^2 * M;

for i = 1:num
    x1 = rand(n,1); % 初值
    for j = 1:i
        c = phi(:,j)' * M * x1;
        x1 = x1 - c * phi(:,j); % 正交化
    end
    x1 = x1 / sqrt( x1' * M * x1 ); % 归一化
    y1 = M * x1;
    k = 1;
    rho1 = x1' * Kbar * x1 / ( x1' * y1 );
    x2b = Kbar \ y1;
% %     
%     x2b = x2b / sqrt(x2b'*M*x2b);
% %     
    y2b = M * x2b;
    rho2 = x2b' * y1 / ( x2b' * y2b );
    y2 = y2b / sqrt( x2b' * y2b );

    while abs(rho2-rho1)/rho2 > tol
        if (k>1e4)
            fprintf('迭代10000次仍未收敛，停止计算\n');
            break;
        end
        k = k + 1;
        y1 = y2;
        rho1 = rho2;
        x2b = Kbar \ y1;
        for j = 1:i
            c = phi(:,j)' * M * x2b;
            x2b = x2b - c * phi(:,j); % 正交化
        end
% %         
%         x2b = x2b / sqrt(x2b'*M*x2b);
% % 
        y2b = M * x2b;
        rho2 = ( x2b' * y1 ) / ( x2b' * y2b );
        y2 = y2b / sqrt( x2b' * y2b );
    end
    lambda(i) = rho2 - alpha^2;
    phi(:,i) = x2b / sqrt( x2b' * y2b );
end

end

function [x,lam,k] = inverse1(K,M,tol)
    n = size(K,1);
    x1 = ones(n,1);
    y1 = M*x1;
    k = 1;
    x2b = K\y1;
    y2b = M*x2b;
    rho1 = x1'*K*x1/(x1'*M*x1);
    rho2 = x2b'*y1/(x1'*y2b);
    y2 = y2b/sqrt( x2b'*y2b );
    while abs(rho2-rho1)/rho2 > tol
        if (k>1e4)
            fprintf('迭代10000次仍未收敛，停止计算\n');
            break;
        end
        k = k + 1;
        y1 = y2;
        x2b = K\y1;
        y2b = M*x2b;
        rho1 = rho2;
        rho2 = x2b'*y1/(x2b'*y2b);
        y2 = y2b/sqrt( x2b'*y2b );
    end
    lam = rho2;
    x = x2b/sqrt(x2b'*y2b);
end