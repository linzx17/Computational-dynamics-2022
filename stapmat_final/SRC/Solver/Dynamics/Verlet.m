% % Verlet法求解时间积分
function Verlet()

global cdata;
global sdata;

K = sdata.SPSTIFF; % 刚度阵
M = sdata.SPMASS; % 质量阵
Calpha = sdata.Calpha;
Cbeta = sdata.Cbeta;
C = Calpha * M + Cbeta * K; % 阻尼阵

fre_max = max(sdata.FRE_EIG);
dtcr = 1 / fre_max / pi; % 临界时间步长
coef1 = 0.8;
dt = coef1 * dtcr;
t_end = cdata.DSTIME;
t = ( 0:dt:t_end )';


% % % 
n_K = size(K,1);
n = length(t);

u = zeros(n_K,n);%位移
v = zeros(n_K,n);%速度
a = zeros(n_K,n);%加速度
FLOAD = sdata.R(:,1); % 节点载荷
FPAR = GetDynamicsPAR(t); % 动载荷系数
Q = FLOAD * FPAR;
a(:,1) = M \ ( Q(:,1) - K * u(:,1) );%初值

for i = 1:n-1
    v2 = v(:,i) + 1/2 * a(:,i) * dt;%v(t+1/2dt)
    u(:,i+1) = u(:,i) + v2*dt;
    a(:,i+1) = M \ ( Q(:,i+1) - K * u(:,i+1) );
    v(:,i+1) = v2 + 1/2 * a(:,i+1) * dt;
end
u_verlet = u';

% % 测试结果
figure
plot(t,u_verlet(:,260),'LineWidth',1.5);
xlabel('time(s)'); ylabel('u_{260}');
title('速度Verlet');
set(gca,'FontSize',16);
% % 

sdata.DYNADT = t;
sdata.DYNADIS = u_verlet;
% % % 

end