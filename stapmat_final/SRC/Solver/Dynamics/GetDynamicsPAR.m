%* *****************************************************************
%* - Function of STAPMAT in Solver/Dynamics phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     计算一个时间序列上的动载荷系数        *
%*                                                                 *
%* - Call procedures:                                              *
%*                          *
%*                                                                 *
%* - Called by :                                                   *
%*                                                        *
%*                                                                 *
%* - Programmed by:                                                *
%*     Zixiong Lin, 2022.05.30                *
%*                                                                 *
%* *****************************************************************
function FPAR = GetDynamicsPAR(t)

global cdata;
global sdata;
DLTYPE = cdata.DLTYPE;
DLPAR = cdata.DLPAR;
% FLOAD = sdata.R(:,1); % 节点载荷

if 2 == DLTYPE
    t1 = DLPAR(1);
    t4 = DLPAR(4);
    w = DLPAR(2);
    phi = DLPAR(3);
    FPAR = sin_dyload(t,t1,t4,w,phi); % 动载荷系数：sin函数
elseif 1 == DLTYPE
    t1 = DLPAR(1);
    t2 = DLPAR(2);
    t3 = DLPAR(3);
    t4 = DLPAR(4);
    FPAR = piecewise_linear_dyload(t,t1,t2,t3,t4); % 动载荷系数：分段线性
end

% DYFLOAD = FPAR * FLOAD; % 节点动载荷向量

end

%% --------------------functions--------------------

% % sin函数式动载荷系数
% % 输入待求时刻t，载荷起止时刻t1,t4，振动角频率w，初始相位phi
% % 输出动载荷系数
function FPAR = sin_dyload(t,t1,t4,w,phi)

n = length(t);
FPAR = zeros(1,n);
for i = 1:n
    if (t(i) >= t1) && (t(i) <= t4)
        FPAR(i) = sin( w * t(i) + phi );
    else
        FPAR(i) = 0;
    end
end

end


% % 分段线性函数式动载荷系数
% % 输入待求时刻t，载荷起止时刻t1,t4，中间时刻t2,t3
% % 输出动载荷系数
function FPAR = piecewise_linear_dyload(t,t1,t2,t3,t4)

n = length(t);
FPAR = zeros(1,n);
for i = 1:n
    if (t(i) < t1) || (t(i) > t4)
        FPAR(i) = 0;
    elseif (t(i) >= t2) && (t(i) <= t3)
        FPAR(i) = 1;
    elseif (t(i) >= t1) && (t(i) < t2)
        FPAR(i) = (t(i)-t1) / (t2-t1);
    elseif (t(i) > t3) && (t(i) <= t4)
        FPAR(i) = (t(i)-t4) / (t3-t4);
    end
end

end