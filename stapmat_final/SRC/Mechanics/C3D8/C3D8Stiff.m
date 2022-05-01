% % Called by GetStiff.m
% % 
% % - Call procedures:
% %        InitC3D8() % 定义每个单元上的节点数和每个节点的自由度
% %        ReadC3D8()
% % 

function C3D8Stiff()

% Init variables of the element
InitC3D8();

% Read Material and Elements
ReadC3D8();

fprintf('Solution phase ...\n\n');

% calculate addresses of diagonal elements
Addres();

% Data check Or Solve
global cdata;
if (cdata.MODEX == 0) 
    cdata.TIM(3,:) = clock;
    cdata.TIM(4,:) = clock;
    cdata.TIM(5,:) = clock;
    return; 
end

% Assemble structure stiffness matrix
Assemble();

end % end of Solid3DStiff


% ----------------------- Functions -----------------------------------

% Init parameters of C3D8 element
function InitC3D8()

global sdata;
sdata.NNODE = 8; % 一个单元上的节点数
sdata.NDOF = 3; % 每个节点的自由度

end % end of function InitC3D8()


% Assemble structure stiffness matrix
function Assemble()

global sdata;
global cdata;
sdata.STIFF = zeros(sdata.NWK, 1, 'double');

NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
E = sdata.E; nu = sdata.nu; NGaussian = sdata.NGaussian; LM = sdata.LM;

for N = 1:NUME
    MTYPE = MATP(N);

    [ng,ksi,eta,zeta,weight] = GetGaussianIntInfo( NGaussian(MTYPE) );

    E0 = E(MTYPE);
    nu0 = nu(MTYPE);
    
end

end % end of function Assemble()

% % 输入数字
% % 输出高斯积分点个数ng，高斯积分点的坐标ksi,eta，权重weight
function [ng,ksi,eta,zeta,weight] = GetGaussianIntInfo(num)
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