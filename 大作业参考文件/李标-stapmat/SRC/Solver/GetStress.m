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

function GetStress(NUM) %NUM是载荷工况

% Different type of element
global cdata;
NUMEG = cdata.NUMEG;

for N = 1:NUMEG  %N是不同的单元组
    NPAR1 = cdata.NPAR(1);
    if (NPAR1 == 1) 
        TrussStress(NUM, N)
    elseif (NPAR1 == 2) 
        PlaneStress(NUM, N)
    elseif (NPAR1 == 3)   %增添了高阶8Q轴对称单元
        AxisymStress(NUM, N)
    else
        error(' *** ERROR *** No Such Element'); 
    end
end



end