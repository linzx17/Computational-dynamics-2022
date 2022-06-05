% % 模态叠加法
function ModalSup()

global cdata;
global sdata;

IOUT = cdata.IOUT;

fprintf('Modal Superposition Starts...\n')
fprintf(IOUT,'Modal Superposition Starts...\n');
time1 = clock;

K = sdata.SPSTIFF; % 刚度阵
M = sdata.SPMASS; % 质量阵
Calpha = sdata.Calpha;
Cbeta = sdata.Cbeta;

w = sdata.FREQUENCY * 2 * pi; % 角频率
phi = sdata.EIGVECTOR;
t_end = cdata.DSTIME; % 求解时间
dt = 2*pi/max(w)/20;
if dt > t_end/10
    dt = t_end/10;
end
t = 0:dt:t_end; % 求解时间序列

fprintf('dt = %e, t_end = %.2f, length(t) = %d\n',dt,t_end,length(t));
fprintf(IOUT,'dt = %e, t_end = %.2f, length(t) = %d\n',dt,t_end,length(t));

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

% % 测试结果
% figure
% output_node = 260;
% output_node = 685;
% plot(t,u_modal(:,output_node),'LineWidth',1.5);
% xlabel('time(s)'); ylabel(['u_{',num2str(output_node),'}']);
% title('模态叠加');
% set(gca,'FontSize',16);
% % 

sdata.DYNADT = t;
sdata.DYNADIS = u_modal;

time2 = clock;
fprintf('Modal Superposition ends. TIME = %.2f\n\n',etime(time2,time1));
fprintf(IOUT,'TIME FOR Modal Superposition  . . . . . . . . . . . . . . = %.2f\n\n',etime(time2,time1));

% % % % % 输出结果，跑完注释掉
% save('E:\lzx2021\dyna\Job-real_large_dyna-1451\real_1451_u_modal.mat','u_modal');
% % % % % 

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

