function Opt_beta()


%������������Ż�
% Get global data
global sdata;
global cdata;
dt = cdata.dt;

K = sdata.STIFF; M = sdata.MASS;
NEQ = sdata.NEQ;
N_w = 10; %�����Ż�Ƶ�ʶ�Ӧ�Ľ���N_w
%beta = sdata.beta;
q = 2*N_w;
if (q > NEQ)
    q = NEQ;
    N_w = NEQ;
end

[w_m,w] = SSI(q,K,M,N_w);%�ӿռ���������ϵͳ�ϸ߽�Ƶ��q��N_w��Ƶ��w

Tmin = 2*pi/w_m;
%�Ż�beta,dt����̫С�����ȡTmin/10
if(dt<Tmin/10)
    c = 0;%�����Ϊ0
    wd = w*(1 - c^2)^0.5;
    ta = (tan(wd*dt))^2;
    B = 1 + 2*c/(w*dt);
    x0 = (2/(B^2))*(1 + B*ta - (1 + 2*B*ta - B*B*ta)^0.5)/(1 + ta);
    beta = 1/x0 - 1/((w*dt)^2) - c/(w*dt)
    cdata.beta = beta;
else
    printf('����ƫ��Ӧ����dt<=%0.4e\n',Tmin/10);
end
    

