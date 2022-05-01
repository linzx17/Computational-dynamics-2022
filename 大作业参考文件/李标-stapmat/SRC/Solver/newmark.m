function newmark(~)
L = 1;
global cdata;
global sdata;
n = cdata.nt;
dt = cdata.dt;
v = 300;%�����ٶ� 300mm/ms���ڳ����п��ƣ���Ϊ���Ƕ����ض��غ������Ƶģ�����ģ�鴦��
gamma = cdata.gamma;
beta = cdata.beta;
NN = sdata.NEQ;
K = sdata.STIFF; M = sdata.MASS; 
C = sparse(NN, NN);
Qq = sdata.R(:,L)/1.2; %�����غɴ�С�������ļ���Ϊÿ��Ԫ0.6�ĺ���
Q = Qq*0;

TIM1 = clock;%��ʱ

%�������
c0 = 1/(beta*dt^2);
c1 = gamma/(beta*dt);
c2 = 1/(beta*dt);
c3 = 1/(2*beta) - 1;
c4 = gamma/beta - 1;
c5 = dt*(gamma/(2*beta) - 1);
c6 = dt*(1 - gamma);
c7 = gamma*dt;

%�γ���Ч�ն���
K_eff = K + c0*M +c1*C;
 K_eff = sparse(K_eff);
%���Ƿֽ�
[L1,D1,P] = ldl(K_eff);

%��ʼ��
xd = zeros(NN,1, 'double');%����������ʹ�õ�λ�ƣ��ٶȣ����ٶȣ��洢����
d1 = xd;
v1 = zeros(NN,1, 'double');
a1 = M\(Q-K*d1); %V0 = 0


%t = zeros(n,1);
%t(1) = 0;



%�������

 for i=1:n
   Q = Qq*sin(i*dt);
    Q_eff = Q + M*(c0*d1 + c2*v1 + c3*a1) + C*(c1*d1 + c4*v1 + c5*a1);
    
    %d2 = K_eff\Q_eff;
    d2 = P*(L1'\(D1\(L1\(P'*Q_eff))));                       
    a2 = c0*(d2 - d1) - c2*v1 - c3*a1;
    v2 = v1 + c6*a1 + c7*a2;
    
    %xd(:,i+1) = d2;
    
    %t(i+1) = t(i) + dt;
    d1 = d2; v1 = v2; a1 = a2;
    
    ddx(i) = d1(1921);
    ddz(i) = d1(1922);
    ttt(i) = i*dt;
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % ���Ϊvtk��ʽ
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (rem(i,1000)==0)%ÿ100�ε�����ɾ����һ�Σ��������ʱ�����ĩֵ
        sdata.DIS(:,1) = d2;
        sdata.V(:,1) = v2;
        GetStress(L);
        timestap = i;
        vtkwrite('newmark',timestap);
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % ���Ϊvtk��ʽ
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
 end
ddxx = (ddx.^2+ddz.^2).^0.5;
figure;
plot(ttt,ddx,'r'),hold on;
plot(ttt,ddz,'k'),hold on;
plot(ttt,ddxx,'m'),hold on;
title('λ��ʱ����Ӧ����')
xlabel('t')
ylabel('x')   

TIM2 = clock;
time(1) = etime(TIM2, TIM1);
fprintf(['\n' ...
    '     �����ͺ���ʱ��  . . . . . . . . . . . . . . = %12.2f\n'] ...
    ,time(1));
end