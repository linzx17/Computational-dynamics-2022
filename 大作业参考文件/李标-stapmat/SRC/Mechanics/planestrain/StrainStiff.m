%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Compute the stiffness matrix of truss                       *
%*                                                                 *
%* - Call procedures:                                              *
%*     TrussStiff.m - InitTruss()                                  *
%*     ./ReadTruss.m - ReadTruss()                                 *
%*     SRC/Mechanics/Addres.m - Addres()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/GetStiff.m                                    *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************

function StrainStiff()

% Init variables of the element
InitStrain();

% Read Material and Elements
planeReadStrain();

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
planeAssemble();




end

% ----------------------- Functions -----------------------------------

% Init parameters of truss element
function InitStrain()
global sdata;
sdata.NNODE = 4;
sdata.NDOF = 2;

end



