% % 输入数字
% % 输出高斯积分点个数ng，高斯积分点的坐标ksi,eta，权重weight
function [ng,ksi,eta,zeta,weight] = GetC3D8GaussianIntInfo(num)
    global sdata;
    if num == 1
        ng = 1;
        ksi = sdata.GC1;
        eta = sdata.GC1;
        zeta = sdata.GC1;
        weight = sdata.GW1;
    elseif num == 2
        ng = 2;
        ksi = sdata.GC2;
        eta = sdata.GC2;
        zeta = sdata.GC2;
        weight = sdata.GW2;
    else
        ng = 3;
        ksi = sdata.GC3;
        eta = sdata.GC3;
        zeta = sdata.GC3;
        weight = sdata.GW3;
    end
end % end of function GetGaussianIntInfo