function nnewmark(~)
L = 1;
global sdata;
global cdata;
v = 300;%激波速度 300mm/ms，在程序中控制，因为这是对于特定载荷情况设计的，不需模块处理
n = cdata.nt;
dt = cdata.dt;
gamma = cdata.gamma;
beta = cdata.beta;
outt = cdata.outt;
NN = sdata.NEQ;
K = sdata.STIFF; M = sdata.MASS;
C = sparse(NN, NN);
Qq = sdata.R(:,L); %调整载荷大小，输入文件中为每单元0.6的合力
Q = Qq;%初始载荷，不同算例手动调整


misesnum = 0;%输出最大米塞斯应力的参数


TIM1 = clock;%计时


q1 = (dt*dt*(1-2*beta)/2);%提前计算系数，减少计算量
q2 = gamma*dt;
q3 = beta*dt*dt;

xd = zeros(NN,1, 'double');%迭代过程中使用的位移，速度，加速度，存储矩阵
xv =zeros(NN,1, 'double');
xa = zeros(NN,1, 'double');

[C1,K1,M1,C2,K2,M2] = split(C,K,M);% 矩阵分解为对角阵和剩余部分之和的形式，M阵为集中质量阵不需要分解
MM1 = diag(M1+q2*C1+q3*K1);
xa(:,1) = M\(Q-K*xd(:,1));%计算初始加速度

for II = 1:n
    
    t = II*dt;
     %Q = Qq*sin(t);%正弦加载的算例的Q更新
    dd0 = xd+dt*xv+q1*xa;
    vv0 = xv+(dt-q2)*xa;
    error11 =10;
    while (1e2>error11)&&(error11>1e-8)
        
        ddn = xd;%记录下来用于判断收敛
        
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
    

    %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    %激波算例更新Q 若不是激波前进算例，全部注释即可,仅对75*300单元适用
    %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
%     CHNOD = sdata.CHNOD;%变化载荷用（包含载荷节点和XY方向坐标，见READFILE）
%     CHNOD = sortrows(CHNOD,3);%按第三行排序，也就是z坐标
%     underpower = (CHNOD(move*2+1,2)-CHNOD(move*2-1,2))*100*2*pi*10/6;%上一单元贡献的力，激波加载使用的变量
%     move = 1;%激波到达单元数
%     if (move<=300)
%         dispower = v*t;
%         if (dispower>CHNOD(move*2+1,2))
%             DD = CHNOD(move*2-1,1);
%             
%             if (DD > 0)
%                 Q(DD) = Qq(DD);%激波经过单元后更新为均匀应力
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
%             xxeq = -1+xxueq/(CHNOD(move*2,2)-CHNOD(move*2-1,2));%激波波头在单元内的等参位置
%             W1 = xxeq^3/6-xxeq^2/4+5/12;
%             W2 = xxeq^3/6+xxeq^2/4-1/12;
%             W3 = xxeq-xxeq^3/3+2/3;
%             power = xxueq*100*2*pi*10;%本单元上合力 *2piR
%             power1 = power*W1/(W1+W2+W3);%各节点力
%             power2 = power*W3/(W1+W2+W3);
%             power3 = power*W2/(W1+W2+W3);
%             DD = CHNOD(move*2-1,1);
%             if (DD > 0)
%                 Q(DD) = power1+Qq(DD)-underpower;%更新各个节点
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
    %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    %更新Q
    %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    
    
    ddx(II) = xd(1921);
    ddz(II) = xd(1922);
    ttt(II) = t;
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % 输出为vtk格式
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (rem(II,outt)==0)%每outt次迭代计算一次应力并输出
        sdata.DIS(:,1) = xd;
        sdata.V(:,1) = xv;
        GetStress(L);
        misesnum = misesnum+1;
        Mises(1,misesnum) = sdata.maxMises(1);%最大Mises应力，用于绘图
        Mises(2,misesnum) = sdata.maxMises(2);
        Mises(3,misesnum) = t;
        vtkwrite('nnewmark',II);
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % 输出为vtk格式
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end

TIM2 = clock;
time(1) = etime(TIM2, TIM1);%显示迭代过程和输出过程的总时间
fprintf(['\n' ...
    '     迭代和后处理时间  . . . . . . . . . . . . . . = %12.2f\n'] ...
    ,time(1));
ddxx = (ddx.^2+ddz.^2).^0.5;
figure;
plot(ttt,ddx,'r'),hold on;
plot(ttt,ddz,'k'),hold on;
plot(ttt,ddxx,'m'),hold on;
title('位移时程响应曲线')
xlabel('t')
ylabel('x')
% figure;
% plot(Mises(3,:),Mises(1,:),'r'),hold on;
% title('最大米塞斯应力')
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
