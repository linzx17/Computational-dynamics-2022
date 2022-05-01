%* *****************************************************************
%* - Function of STAPMAT in Solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To solve finite element static equilibrium equations        *
%*                                                                 *
%* - Call procedures:                                              *
%*     ./LDLTFactor.m            - LDLTFactor()                    *
%*     Solve.m                   - Stiff2Sparse()                  *
%*     ./ColSol.m                - ColSol()                        *  
%*     Solve.m                   - WriteDis()                      *
%*     SRC/Mechanics/GetStress.m - GetStress()                     *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function Solve()

global cdata;
global sdata;

NEQ = sdata.NEQ;
NLCASE = cdata.NLCASE;
MODEX = cdata.MODEX;
NPAR1 = cdata.NPAR(1);
betamode = cdata.betamode;
sdata.DIS = zeros(NEQ, NLCASE, 'double');
sdata.V = zeros(NEQ, NLCASE, 'double');
sdata.STRAIN = zeros(NEQ, NLCASE, 'double');
sdata.STRESS = zeros(NEQ, NLCASE, 'double');

% The pre-process of Solution
% MODEX = 1, 稀疏求解     
% MODEX = 2, nnewmark

if (NPAR1 == 3)
    change();%针对轴对称单元的载荷处理
end
if (MODEX > 1)
    if (betamode == 1)
        Opt_beta();
    end
end
for L = 1:NLCASE
    if (MODEX == 1)
        Qq = sdata.R(:,L);
        sdata.DIS(:,L) = sdata.STIFF \ Qq;
        %   Print displacements
        WriteDis(L);
    
    %   Calculation of stresses
        GetStress(L);
        
        vtkwrite('PARA',0);
    elseif (MODEX == 2)
       
        nnewmark(L);
        
    elseif (MODEX == 3)

        newmark(L);
        
    end
    
end


end

% ----------------------- Functions -----------------------------------


% Print Displacements
function WriteDis(NUM)

% Get global data
global cdata;
global sdata;
IOUT = cdata.IOUT;
NUMNP = cdata.NUMNP;
DIS = sdata.DIS(:, NUM); ID = sdata.ID;

fprintf(IOUT, '\n\n LOAD CASE %3d', NUM);
fprintf(IOUT, ['\n\n D I S P L A C E M E N T S\n' ...
    '\n       NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT\n']);

D = zeros(3, 1, 'double');
for II = 1:NUMNP
    D(:) = 0;
    if (ID(1, II) ~= 0) D(1) = DIS(ID(1, II)); end
    if (ID(2, II) ~= 0) D(2) = DIS(ID(2, II)); end
    if (ID(3, II) ~= 0) D(3) = DIS(ID(3, II)); end
    
    fprintf(IOUT, ' %10d        %18.6e%18.6e%18.6e\n', II, D(1), D(2), D(3));
end

end

function change()
global sdata;
global cdata;
NLOAD = cdata.NLOAD;
CHNOD = sdata.CHNOD;
Q = sdata.R(:,1);
for L = 1:NLOAD
    II = CHNOD(L,1);
    if (II > 0)
        Q(II) = Q(II)*2*pi*CHNOD(L,2);
    end
end
sdata.R(:,1) = Q;
end