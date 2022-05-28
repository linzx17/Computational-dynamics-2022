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
IOUT = cdata.IOUT;

NUM_EIG = 10; % 求解前NUM_GIE阶的特征值
if NUM_EIG > NEQ
    NUM_EIG = NEQ;
end

% sdata.EIGVALUE = zeros(NUM_EIG, 1, 'double');
% sdata.EIGVECTOR = zeros(NEQ, NUM_EIG, 'double');

SPSTIFF = Matrix2Sparse(sdata.STIFF);
SPMASS = Matrix2Sparse(sdata.MASS);
SPMASS_full = full(SPMASS);%去掉约束的刚度阵，总体刚度阵
SPSTIFF_full = full(SPSTIFF);

% % % 使用eig检查结果，可注释
[V,D] = eig(SPSTIFF_full,SPMASS_full);
D_sqrt = sqrt(diag(D))./2.0./pi;
fprintf(IOUT, 'Eigenvalues computed using eig:\n');
for iout = 1:NUM_EIG
    fprintf(IOUT, '%e     ',D(iout,iout));
end
fprintf(IOUT, '\n\n');
fprintf(IOUT, 'Frequency computed using eig:\n');
for iout = 1:NUM_EIG
    fprintf(IOUT, '%e     ',D_sqrt(iout));
end
fprintf(IOUT, '\n\n');
% % % 

tol = 1e-6; % 收敛判据
% [x,lam,k] = inverse1(full(SPSTIFF),full(SPMASS),tol);
% f_1 = sqrt(lam)./2.0./pi;%[Hz]
[lambda,phi] = inverse(full(SPSTIFF),full(SPMASS),NUM_EIG,tol); % 特征值lambda和特征向量phi
f_all = sqrt(lambda)./2.0./pi;%[Hz]

sdata.EIGVALUE = lambda;
sdata.EIGVECTOR = phi;
sdata.FREQUENCY = f_all;

fprintf(IOUT, 'Eigenvalues computed using stapmat:\n');
for iout = 1:NUM_EIG
    fprintf(IOUT, '%e     ',lambda(iout));
end
fprintf(IOUT, '\n\n');
fprintf(IOUT, 'Frequency computed using stapmat:\n');
for iout = 1:NUM_EIG
    fprintf(IOUT, '%e     ',f_all(iout));
end
fprintf(IOUT, '\n\n');
fprintf(IOUT, 'Eigen vector computed using stapmat:\n');
for jout = 1:size(phi,1)
    for iout = 1:NUM_EIG
        fprintf(IOUT, '%e     ',phi(jout,iout));
    end
    fprintf(IOUT, '\n');
end
fprintf(IOUT, '\n\n');

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

% if abs(det(K)) < 1e-6
if min( abs( eig(K) ) ) < 1e-6
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