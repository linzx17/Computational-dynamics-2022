%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     Calculation of strain and stress                            *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Mechanics/Truss/TrussStress.m   - TrussStress()         *       
%*                                                                 *
%* - Called by :                                                   *
%*     ./Solver.m                                                  *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function GetStress(NUM) % NUM:载荷工况数

% Different type of element
global cdata;
NUMEG = cdata.NUMEG;
IOUT = cdata.IOUT;

for N = 1:NUMEG % 单元组数
    NPAR1 = cdata.NPAR(1);
    if (NPAR1 == 1) 
        TrussStress(NUM, N);
    elseif (NPAR1 == 2)
        PlaneStress(NUM, N);
    elseif (NPAR1 == 3)
        C3D8Stress(NUM, N);
    elseif (NPAR1 == 4)
        C3D20Stress(NUM, N);
    else 
        error(' *** ERROR *** No Such Element');
    end 
end



end