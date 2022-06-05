% % 广义alpha法求时间积分
function GenAlpha()

global cdata;
global sdata;

IOUT = cdata.IOUT;

fprintf('Generalized Alpha Starts...\n')
fprintf(IOUT,'Generalized Alpha Starts...\n');
time1 = clock;

K = sdata.SPSTIFF; % 刚度阵
M = sdata.SPMASS; % 质量阵

Calpha = sdata.Calpha;
Cbeta = sdata.Cbeta;
C = Calpha * M + Cbeta * K; % 阻尼阵

% fre_max = max(sdata.FRE_EIG);
fre_max = max(sdata.FREQUENCY);
dtcr = 1 / fre_max / pi; % 临界时间步长
coef1 = 0.8;
dt = coef1 * dtcr;
t_end = cdata.DSTIME;
t = ( 0:dt:t_end )';
rho = 1;

fprintf('dt = %e, t_end = %.2f, length(t) = %d\n',dt,t_end,length(t));
fprintf(IOUT,'dt = %e, t_end = %.2f, length(t) = %d\n',dt,t_end,length(t));

Q = sdata.R(:,1); % 节点载荷

u_alpha = generalized_alpha(K,M,C,Q,t,rho);

% % 测试结果
% figure
% % output_node = 260;
% output_node = 685;
% plot(t,u_alpha(:,output_node),'LineWidth',1.5);
% xlabel('time(s)'); ylabel(['u_{',num2str(output_node),'}']);
% title('广义alpha');
% set(gca,'FontSize',16);
% % 

sdata.DYNADT = t;
sdata.DYNADIS = u_alpha;

time2 = clock;
fprintf('Generalized Alpha ends. TIME = %.2f\n\n',etime(time2,time1));
fprintf(IOUT,'TIME FOR Generalized Alpha  . . . . . . . . . . . . . . = %.2f\n\n',etime(time2,time1));

% % % % % 输出结果，跑完注释掉
% save('E:\lzx2021\dyna\Job-real_large_dyna-1451\real_1451_u_alpha.mat','u_alpha');
% % % % % 

end

function u_alpha = generalized_alpha(K,M,C,Q0,ta,rho)

n_K = size(K,1);
nta = length(ta);
dt = ta(2) - ta(1);

alphaf = rho/(rho+1);
alpham = (2*rho-1)/(rho+1);
beta = (1-alpham+alphaf)^2/4;
gamma = 1/2 - alpham + alphaf;
ck = 1-alphaf;
c0 = (1-alpham)/(beta*dt^2);
c1 = ck*gamma/(beta*dt);
c2 = dt*c0;
c3 = c2*dt/2 - 1;
c4 = ck*gamma/beta - 1;
c5 = ck*(gamma/2*beta-1)*dt;

c11 = 1/(beta*dt^2);
c12 = 1/(beta*dt);
c13 = 1/(2*beta)-1;
c21 = gamma/(beta*dt);
c22 = 1-gamma/beta;
c23 = (1-gamma/(2*beta))*dt;

ua = zeros(n_K,nta);
va = zeros(n_K,nta);
aa = zeros(n_K,nta);
FPAR1 = GetDynamicsPAR(ta(1)); % 动载荷系数
Q = FPAR1 * Q0;
aa(:,1) = M\(Q-K*ua(:,1));
Kbar = ck*K + c0*M + c1*C;
for i = 1:nta-1
    FPAR = GetDynamicsPAR(ta(i+1));
    Q = FPAR * Q0;
    Qbar = Q - alphaf*K*ua(:,i) + M*(c0*ua(:,i)+c2*va(:,i)+c3*aa(:,i)) + C*(c1*ua(:,i)+c4*va(:,i)+c5*aa(:,i));
    ua(:,i+1) = Kbar\Qbar;
    aa(:,i+1) = c11*(ua(:,i+1)-ua(:,i)) - c12*va(:,i) - c13*aa(:,i);
    va(:,i+1) = c21*(ua(:,i+1)-ua(:,i)) + c22*va(:,i) + c23*aa(:,i);
end

u_alpha = ua';

end