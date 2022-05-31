%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Forming the stiffness matrix                                *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Mechanics/Truss/TrussStiff.m - TrussStiff()             *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************

function GetStiff()
% Get global variables
global cdata;
global sdata;

% Read the type of element
IIN = cdata.IIN;
IOUT = cdata.IOUT;
IDAT_ANIM = cdata.IDAT_ANIM;
fprintf(IOUT, '\n\n E L E M E N T   G R O U P   D A T A\n');

NUMEM = int64(0);%初始化单元总数
sdata.NUMEGEM = zeros(cdata.NUMEG,1,'int8');%初始化各单元组中单元序号
sdata.ELEII = zeros(cdata.NUMEM, sdata.NNODE, 'double'); % 单元内各个节点编号信息

for N = 1:cdata.NUMEG % 遍历单元组
    if (N ~= 1)
        fprintf(IOUT,'\n');
    end
    tmp = str2num(fgetl(IIN));
    for I = 1:length(tmp)
        cdata.NPAR(I) = tmp(I);
        % % NPAR(1) 单元类型：1.桁架 2.平面应力 3.三维实体C3D8 4.三维实体C3D20
        % % NPAR(2) 该单元组内的单元数
        % % NPAR(3) 材料/截面属性数
        % % NPAR(4) 该单元类型采用的高斯积分阶数
    end
    if cdata.NPAR(1) == 4
        cdata.NPAR(4) = 3; % 如果单元类型为C3D20，则高斯积分阶数改为3，如果小于3，则无法正确应力磨平
    end
    
    sdata.NUMEGEM(N,1) = NUMEM + cdata.NPAR(2);

    NUMEM = NUMEM + cdata.NPAR(2);    

    fprintf(IOUT, '\n\n E L E M E N T   D E F I N I T I O N\n');
    fprintf(IOUT, ['\n ELEMENT TYPE  . . . . . . . . . . . . .( NPAR(1) ) . . = %10d\n' ...
        '     EQ.1, TRUSS ELEMENTS\n' ...
        '     EQ.2, PLANE STRESS\n' ...
        '     EQ.3, SOLID C3D8\n' ...
        '     EQ.4, SOLID C3D20\n' ...
        ' NUMBER OF ELEMENTS. . . . . . . . . . .( NPAR(2) ) . . = %10d\n'], ...
        cdata.NPAR(1), cdata.NPAR(2));

%   Different kinds of element
    NPAR1 = cdata.NPAR(1);
    if (NPAR1 == 1) % 桁架单元
        TrussStiff();
    elseif (NPAR1 == 2) % 平面应力单元
        PlaneStiff();
    elseif (NPAR1 == 3) % 三维实体单元C3D8
        C3D8Stiff(N);
    elseif (NPAR1 == 4) % 三维实体单元C3D20
        C3D20Stiff(N);
    else 
        error(' *** ERROR *** No Such Element');
    end
    
end

if (cdata.NUMEM ~= NUMEM)
    fprintf(IOUT,'\nTotal Number of elements is error!');
    return
end

fprintf(IDAT_ANIM,'ZONE F=FEPOINT N = %5d E = %5d ET=BRICK C=CYAN \n',cdata.NUMNP,cdata.NUMEM);

sdata.SPSTIFF = Matrix2Sparse(sdata.STIFF);
sdata.SPMASS = Matrix2Sparse(sdata.MASS);

% % % 测试：质量阵使用直接组装的协调质量阵
% sdata.SPMASS = sdata.MASS_DIRCTION;
% % % 

end

% % ----------------------- Functions -----------------------------------

% Convert the stiff vector to a sparse stiff matrix
function SPMatrix = Matrix2Sparse(A)

global sdata;
% A = sdata.STIFF;
MAXA = sdata.MAXA;
NEQ = sdata.NEQ;
NWK = sdata.NWK;
IIndex = zeros(NWK*2-NEQ, 1);
JIndex = IIndex;
STIFF = IIndex;

NUM = 1;
NUMC = 0;
for N = 1:NEQ
    KU = MAXA(N + 1) - MAXA(N);
    for L = 1:KU
        IIndex(NUM) = N;
        JIndex(NUM) = N - L + 1;
        STIFF(NUM) = A(NUM);
        NUM = NUM + 1;
        if (L == 1)
            NUMC = NUMC + 1;
            continue;
        end
        SYMN = NUM-1 - NUMC + NWK;
        IIndex(SYMN) = N - L + 1;
        JIndex(SYMN) = N;
        STIFF(SYMN) = A(NUM-1);
    end
end

SPMatrix = sparse(IIndex, JIndex, STIFF, NEQ, NEQ);
end
