% % 模态叠加法
function ModalSup()

global cdata;
global sdata;

K = sdata.SPSTIFF; % 刚度阵
M = sdata.SPMASS; % 质量阵
Calpha = sdata.Calpha;
Cbeta = sdata.Cbeta;

w = sdata.FREQUENCY * 2 * pi; % 角频率
phi = sdata.EIGVECTOR;
t_end = cdata.DSTIME; % 求解时间
dt = 2*pi/max(w)/20;
t = 0:dt:t_end; % 求解时间序列

n_phi = size(phi,2);
Q0 = sdata.R(:,1); % 节点载荷
% Q = GetDynamicsLoad(t(ti));
phi_Q0 = phi' * Q0; % 阵型坐标下的载荷
% FPAR = GetDynamicsPAR(t); % 时间序列上的动载荷系数
% phi_Q = phi_Q0 * FPAR; % 动载荷
phi_Q = phi_Q0;
for i = 1:n_phi
    [time,q] = ode45(@(t,y)modelfunction(t, y,w(i), phi_Q(i,:), Calpha, Cbeta ), t, [0 0]);
    finalq(:,i) = q(:,1);
end


u_modal = phi * finalq';
u_modal = u_modal';

% % % 测试结果
% figure
% plot(t,u_modal(:,260),'LineWidth',1.5);
% xlabel('time(s)'); ylabel('u_{260}');
% title('模态叠加');
% set(gca,'FontSize',16);
% % % 

sdata.DYNADT = t;
sdata.DYNADIS = u_modal;

end

%%微分方程
function x = modelfunction(t,y,w,F,alpha,beta)
x = zeros(2,1);
x(1) = y(2);
% x(2) = -0.001*w^2*y(2)-w^2*y(1)+F*sin(6*pi*t);
% x(2) = -w^2*y(1)+F;
FPAR = GetDynamicsPAR(t);
x(2) = -(alpha + beta*w^2) * y(2) - w^2*y(1) + F * FPAR;
end

