%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     To assemble element stiffness into global stiffness         *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     ./Truss/ReadTruss.m - Assemble()                            *
%*                                                                 *
%* - Programmed in Fortran 90 by Xiong Zhang                       *
%*                                                                 *
%* - Adapted to Matlab by:                                         *
%*     LeiYang Zhao, Yan Liu, Computational Dynamics Group,        *
%*     School of Aerospace Engineering, Tsinghua University,       *
%*     2019.02.22                                                  *
%*                                                                 *
%* *****************************************************************

function ADDBAN(S, LM,N)

% Get global data
global sdata;
STIFF = sdata.STIFF;NUME = sdata.NUME;
ND = sdata.NDOF * sdata.NNODE;
N1 =ND*ND*(N-1);
for J = 1:ND
    JJ = LM(J);
    for I = 1:ND
        II = LM(I);
        STIFF(N1+ND*(J-1)+I,1) =  II;
        STIFF(N1+ND*(J-1)+I,2) =  JJ;
        STIFF(N1+ND*(J-1)+I,3) =  S(I, J);
    end
end
if (N==NUME)
    STIFF(STIFF(:,1)==0,:)=[];
    STIFF(STIFF(:,2)==0,:)=[];
    STIFF = sparse(STIFF(:,1),STIFF(:,2),STIFF(:,3),sdata.NEQ,sdata.NEQ);
    sdata.STIFF = STIFF;
else
    sdata.STIFF = STIFF;
end

end