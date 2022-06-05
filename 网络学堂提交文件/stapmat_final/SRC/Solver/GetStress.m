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
global sdata;
NUMEG = cdata.NUMEG;
IOUT = cdata.IOUT;
IDAT_CURV = cdata.IDAT_CURV;

sdata.NODE_FLAG = zeros(cdata.NUMNP,cdata.NLCASE,'double');
sdata.STRAIN = zeros(cdata.NUMNP,6,cdata.NLCASE,'double');
sdata.STRESS = zeros(cdata.NUMNP,6,cdata.NLCASE,'double');

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

%应力磨平
for j = 1:cdata.NLCASE
    for i = 1:cdata.NUMNP
        if sdata.NODE_FLAG(i) ~= 1
            sdata.STRAIN(i,1:6,j) = sdata.STRAIN(i,1:6,j)./sdata.NODE_FLAG(i);
            sdata.STRESS(i,1:6,j) = sdata.STRESS(i,1:6,j)./sdata.NODE_FLAG(i);
        end
        fprintf(IDAT_CURV,'%15d%',i);
        for ii = 1:3
            fprintf(IDAT_CURV,'%15.6E%',sdata.DISP(i,ii,j));
        end
        for ii = 1:6
            fprintf(IDAT_CURV,'%15.6E%',sdata.STRAIN(i,ii,j));%应变
        end
        for ii = 1:6
            fprintf(IDAT_CURV,'%15.6E%',sdata.STRESS(i,ii,j));%应力
        end
        fprintf(IDAT_CURV,'\n');
    end
    fprintf(IDAT_CURV,'\n');
end


end