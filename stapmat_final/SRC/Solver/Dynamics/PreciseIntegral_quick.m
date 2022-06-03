% % 具有矩阵稀疏化的精细积分法
% % 第二版精细积分
function PreciseIntegral_quick()

global cdata;
global sdata;
IOUT = cdata.IOUT;
fprintf('Precise Intgral Sparse Start...\n')
fprintf(IOUT,'Precise Intgral Sparse Start...\n');
time1 = clock;

K = full(sdata.SPSTIFF); % 刚度阵
M = full(sdata.SPMASS); % 质量阵
Calpha = sdata.Calpha;
Cbeta = sdata.Cbeta;
C = Calpha * M + Cbeta * K; % 阻尼阵

% dert_t = 1 / max(sdata.FREQUENCY) / 20;%时间域离散步长
% dert_t = 0.05;
dert_t = 0.005;
tn = cdata.DSTIME;%求解的终止时间
t = 0:dert_t:tn;%时间域离散节点向量,必须注意，求解时间必须从0时刻开始否则会有bug
t_stepnum = length(t);%时间离散区间数目

fprintf('dt = %e, t_end = %.2f, length(t) = %d\n\n',dert_t,tn,length(t));
fprintf(IOUT,'dt = %e, t_end = %.2f, length(t) = %d\n\n',dert_t,tn,length(t));

a_dimension = size(M,2);%空间离散变量a的维数
a = zeros(a_dimension,t_stepnum);%用来存放求解的每个时刻a的结果,同一时刻的结构在同一列中
p = zeros(a_dimension,t_stepnum);%增维向量,同一时刻的结构在同一列中
N = 20;%矩阵指数求解时间区间的细分阶数，即将时间区间细分为2^N份

%% 初值写入
a0 = zeros(a_dimension,1);%a的初值
a_dt0 = zeros(a_dimension,1);%a对时间一阶导数的初值
a(:,1) = a0;%a初值的写入；
p(:,1) = M*a_dt0 + C*a0./2;%p初值写入
x = [a;p];%求解的增维变量,同一时刻的结构在同一列中

%% 求解方程构建,R矩阵有专门函数产生
%注意检查质量阵的奇异性，如果是奇异矩阵则求解会出错
I = eye(size(M,2));
M_ni = I/M;
A1 = -M_ni * C ./2;
A2 = C * M_ni * C./4 - K;
A3 = -C * M_ni./2;
A4 = M_ni;
H = [A1 A4; A2  A3];
I_H = eye(size(H,1));
H_ni = H\I_H;
%% 方程求解
T = MatrixIndex(H,dert_t,N);
for i =1:1:t_stepnum-1
%     tic;
%     x(:,i+1) = T * x(:,i) + integralGass(t(i),t(i+1),H);
x(:,i+1) =  T*(x(:,i) + H_ni * ( R(t(i)) + H_ni * R(t(i+1)) ) ) - H_ni * ( R(t(i)) + H_ni * R(t(i+1)) + R(t(i+1)).*dert_t);
%     toc
end

%% 结果回代
a = x(1:a_dimension,:);
p = x(a_dimension+1:end,:);

u_precise_quick = a';

sdata.DYNADT = t;
sdata.DYNADIS = a';

% % 测试结果
% figure
% % output_node = 260;
% output_node = 685;
% plot(t,a(output_node,:),'LineWidth',1.5);
% xlabel('time(s)'); ylabel(['u_{',num2str(output_node),'}']);
% title('精细积分');
% set(gca,'FontSize',16);
% % 

time2 = clock;
fprintf('Precise Intgral Sparse ends. TIME = %.2f\n\n',etime(time2,time1));
fprintf(IOUT,'TIME FOR Precise Intgral Sparse  . . . . . . . . . . . . . . = %.2f\n\n',etime(time2,time1));

% % % % % 输出结果，跑完注释掉
% save('E:\lzx2021\dyna\Job-real_large_dyna-1451\real_1451_u_precise_quick.mat','u_precise_quick');
% % % % % 

end

%% --------------------functions--------------------

function [T] = MatrixIndex(H,dert_t,N)
%求解矩阵指数的函数，H为系数矩阵，dert_t为时间区间，N为细分时间区间的二次方指数，即区间份数m = 2^N
%这里Ta的求解使用泰勒展开，保留四项
  kapa = 1E-25;%稀疏化指标
  dert_tao = dert_t / 2^N;
  
  B =  H .*dert_tao;  
  I = eye(size(B,1));
  Ta = B + B * B./2 *(I + B./3 + B * B./12) ;

  Ta = clear_zero(full(Ta),kapa);
  for i = 1:1:N
      Ta = 2.* Ta + Ta * Ta;
      Ta = clear_zero(Ta,kapa);
  end
  T = I + Ta;
end

function [r] = R(t)
%载荷函数，给时间t，为一个标量,返回一个2n×1的列向量，物理载荷Q在这个函数里面定义
%Q = [-sin(t); 0.5*sin(t)];%载荷矩阵同一时刻的载荷在同一列中
%   global R_total;
%   r = R_total;%zeros(1800,1800);%[zeros(length(Q),1); Q];
global sdata;
Q0 = sdata.R(:,1); % 节点载荷
FPAR = GetDynamicsPAR(t);
Q = Q0 * FPAR;
r = [zeros(size(Q));Q]; % 增广载荷列向量
end

function [result] = integralGass(tk,tk1,H)
%精细积分法，积分部分的求解函数，使用高斯积分，tk为起始时刻，tk1为终止时刻,n为高斯点个数，采用三节点高斯积分
  A = [8/9  5/9  5/9];%积分系数
  y = [0  -0.6^0.5  0.6^0.5];
  n = size(H,2)/2;%H矩阵的维数的一半，即原问题增维之前的维数
  
  N = 20;%矩阵指数区间细分指数，区间数为2^N
  dert_t = tk1 - tk;%求解的时间区间
  result = 0;
  for i = 1:1:length(A)
      r = R(tk + 0.5 * dert_t *(1 + y(i)));%积分节点上对应的载荷向量
      flag = find_unzero(r);
      if flag == 2*n + 1
        continue    
      else
        Ti = MatrixIndex(H,0.5*dert_t*(1 - y(i)),N);
        Ti_low = Ti(:, flag :2*n);%用于加速计算的降低维数的Ti
        R_low = r(flag:2*n, :);%用于加速计算的降低维数的Ri
        result = result + 0.5 * dert_t * A(i) .* Ti_low * R_low; 
      end
  end   
end

function [flag] = find_unzero(R)
%寻找载荷R向量第一个非零元的函数，返回值flag为第一个非零元的行号，如果全为零则返回2n+1
  n = length(R)/2;
  R_find = R(n+1:end,:);%用来找第一个非零元的载荷向量子块，因为R的前N项都是零元，所以直接从N +1开始找
  unzero_num = find(R_find~=0);%非零元的序号，是一个列向量
  if length(unzero_num) == 0
      flag = 2*n + 1;
  else
      flag = n + unzero_num(1);%载荷项第一个非零行的位置
  end
end

function [Ta] = clear_zero(Ta,kapa)
%用于矩阵指数迭代过程中误差项的归零，kapa是归零指标，一般设置为1E-25;
 n = size(Ta,2)/2;
 Ta11 = Ta(1:n,1:n);
 Ta12 = Ta(1:n,n+1:2*n);
 Ta21 = Ta(n+1:2*n,1:n);
 Ta22 = Ta(n+1:2*n,n+1:2*n);
 
 abs_Ta11 = abs(Ta11);
 abs_Ta12 = abs(Ta12);
 abs_Ta21 = abs(Ta21);
 abs_Ta22 = abs(Ta22);
 
 max_num11 = max(max(abs_Ta11));
 max_num12 = max(max(abs_Ta12));
 max_num21 = max(max(abs_Ta21));
 max_num22 = max(max(abs_Ta22));
 
 Ta11(abs_Ta11 <= max_num11*kapa) = 0;
 Ta12(abs_Ta12 <= max_num12*kapa) = 0;
 Ta21(abs_Ta21 <= max_num21*kapa) = 0;
 Ta22(abs_Ta22 <= max_num22*kapa) = 0;
 
 Ta = [Ta11 Ta12; Ta21 Ta22;]; 
end