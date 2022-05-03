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

% Read the type of element
IIN = cdata.IIN;
IOUT = cdata.IOUT;
fprintf(IOUT, '\n\n E L E M E N T   G R O U P   D A T A\n');

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
        C3D8Stiff();
    elseif (NPAR1 == 4) % 三维实体单元C3D20
        C3D20Stiff();
    else 
        error(' *** ERROR *** No Such Element');
    end
    
end

end
