function nnewmark(~)
L = 1;
global sdata;
global cdata;
v = 300;%�����ٶ� 300mm/ms���ڳ����п��ƣ���Ϊ���Ƕ����ض��غ������Ƶģ�����ģ�鴦��
n = cdata.nt;
dt = cdata.dt;
gamma = cdata.gamma;
beta = cdata.beta;
outt = cdata.outt;
NN = sdata.NEQ;
K = sdata.STIFF; M = sdata.MASS;
C = sparse(NN, NN);
Qq = sdata.R(:,L); %�����غɴ�С�������ļ���Ϊÿ��Ԫ0.6�ĺ���
Q = Qq;%��ʼ�غɣ���ͬ�����ֶ�����


misesnum = 0;%����������˹Ӧ���Ĳ���


TIM1 = clock;%��ʱ


q1 = (dt*dt*(1-2*beta)/2);%��ǰ����ϵ�������ټ�����
q2 = gamma*dt;
q3 = beta*dt*dt;

xd = zeros(NN,1, 'double');%����������ʹ�õ�λ�ƣ��ٶȣ����ٶȣ��洢����
xv =zeros(NN,1, 'double');
xa = zeros(NN,1, 'double');

[C1,K1,M1,C2,K2,M2] = split(C,K,M);% ����ֽ�Ϊ�Խ����ʣ�ಿ��֮�͵���ʽ��M��Ϊ������������Ҫ�ֽ�
MM1 = diag(M1+q2*C1+q3*K1);
xa(:,1) = M\(Q-K*xd(:,1));%�����ʼ���ٶ�

for II = 1:n
    
    t = II*dt;
     %Q = Qq*sin(t);%���Ҽ��ص�������Q����
    dd0 = xd+dt*xv+q1*xa;
    vv0 = xv+(dt-q2)*xa;
    error11 =10;
    while (1e2>error11)&&(error11>1e-8)
        
        ddn = xd;%��¼���������ж�����
        
        Q1 = Q+M2*xa+C2*xv+K2*xd-K1*dd0-C1*vv0;
        
        xa = Q1./MM1;
        xd = dd0+q3*xa;
        xv = vv0+q2*xa;
        
        
        error11 = abs(norm(ddn,1)-norm(xd,1));
        if (error11 ~= 0)
            error11 = error11/abs(norm(xd,1));
        end
        
    end
    
    if (error11>1)||isnan(error11)
        error(' *** ERROR *** boom');
    end
    

    %��������������������������������������������������������������
    %������������Q �����Ǽ���ǰ��������ȫ��ע�ͼ���,����75*300��Ԫ����
    %��������������������������������������������������������������
%     CHNOD = sdata.CHNOD;%�仯�غ��ã������غɽڵ��XY�������꣬��READFILE��
%     CHNOD = sortrows(CHNOD,3);%������������Ҳ����z����
%     underpower = (CHNOD(move*2+1,2)-CHNOD(move*2-1,2))*100*2*pi*10/6;%��һ��Ԫ���׵�������������ʹ�õı���
%     move = 1;%�������ﵥԪ��
%     if (move<=300)
%         dispower = v*t;
%         if (dispower>CHNOD(move*2+1,2))
%             DD = CHNOD(move*2-1,1);
%             
%             if (DD > 0)
%                 Q(DD) = Qq(DD);%����������Ԫ�����Ϊ����Ӧ��
%             end
%             DD = CHNOD(move*2,1);
%             if (DD > 0)
%                 Q(DD) = Qq(DD);
%             end
%             DD = CHNOD(move*2+1,1);
%             if (DD > 0)
%                 Q(DD) = Qq(DD);
%             end
%             move = move+1;
%         end
%         if (move<=300)
%             xxueq = dispower-CHNOD(move*2-1,2);
%             xxeq = -1+xxueq/(CHNOD(move*2,2)-CHNOD(move*2-1,2));%������ͷ�ڵ�Ԫ�ڵĵȲ�λ��
%             W1 = xxeq^3/6-xxeq^2/4+5/12;
%             W2 = xxeq^3/6+xxeq^2/4-1/12;
%             W3 = xxeq-xxeq^3/3+2/3;
%             power = xxueq*100*2*pi*10;%����Ԫ�Ϻ��� *2piR
%             power1 = power*W1/(W1+W2+W3);%���ڵ���
%             power2 = power*W3/(W1+W2+W3);
%             power3 = power*W2/(W1+W2+W3);
%             DD = CHNOD(move*2-1,1);
%             if (DD > 0)
%                 Q(DD) = power1+Qq(DD)-underpower;%���¸����ڵ�
%             end
%             DD = CHNOD(move*2,1);
%             if (DD > 0)
%                 Q(DD) = power2;
%             end
%             DD = CHNOD(move*2+1,1);
%             if (DD > 0)
%                 Q(DD) = power3;
%             end
%         end
%     end
    %��������������������������������������������������������������
    %����Q
    %��������������������������������������������������������������
    
    
    ddx(II) = xd(1921);
    ddz(II) = xd(1922);
    ttt(II) = t;
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % ���Ϊvtk��ʽ
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (rem(II,outt)==0)%ÿoutt�ε�������һ��Ӧ�������
        sdata.DIS(:,1) = xd;
        sdata.V(:,1) = xv;
        GetStress(L);
        misesnum = misesnum+1;
        Mises(1,misesnum) = sdata.maxMises(1);%���MisesӦ�������ڻ�ͼ
        Mises(2,misesnum) = sdata.maxMises(2);
        Mises(3,misesnum) = t;
        vtkwrite('nnewmark',II);
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % ���Ϊvtk��ʽ
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end

TIM2 = clock;
time(1) = etime(TIM2, TIM1);%��ʾ�������̺�������̵���ʱ��
fprintf(['\n' ...
    '     �����ͺ���ʱ��  . . . . . . . . . . . . . . = %12.2f\n'] ...
    ,time(1));
ddxx = (ddx.^2+ddz.^2).^0.5;
figure;
plot(ttt,ddx,'r'),hold on;
plot(ttt,ddz,'k'),hold on;
plot(ttt,ddxx,'m'),hold on;
title('λ��ʱ����Ӧ����')
xlabel('t')
ylabel('x')
% figure;
% plot(Mises(3,:),Mises(1,:),'r'),hold on;
% title('�������˹Ӧ��')
% xlabel('t')
% ylabel('x')

end



function [C1,K1,M1,C2,K2,M2] = split(C,K,M)
C1 = diag(diag(C));
K1 = diag(diag(K));
M1 = diag(diag(M));
C2 = C1-C;
K2 = K1-K;
M2 = M1-M;
end
