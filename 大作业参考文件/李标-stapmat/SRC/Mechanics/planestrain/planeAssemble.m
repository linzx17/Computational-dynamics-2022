%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Assemble structure stiffness matrix                         *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Mechanics/ADDBAN.m - ADDBAN()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/StrainStiff.m                                 *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************
function planeAssemble()
global sdata;
global cdata;

sdata.STIFF = zeros(sdata.NWK, 1, 'double');
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
EN = sdata.E; AREA = sdata.AREA; LM = sdata.LM; NINT = sdata.NINT;
PRN = sdata.PR; IOUT = cdata.IOUT;

for N = 1:NUME
    S = zeros(8, 8, 'double');
    XX = zeros(2, 4, 'double');
    B = zeros(3, 8, 'double'); 
    DB = zeros(3, 'double');
    XG = zeros(4, 4, 'double');
    WGT = zeros(4, 4, 'double');
    MTYPE = MATP(N);
    NU = PRN(MTYPE);
    E = EN(MTYPE);
    for i=1:2
        for j=1:4
            XX(i,j) = XYZ((j-1)*2+i,N);
        end
    end
    
if N==4 %¶Ïµã
end
    
    XG=[0.0  -0.5773502691896  -0.7745966692415  -0.8611363115941;...
        0.0   0.5773502691896   0.0              -0.3399810435849;...
        0.0   0.0               0.7745966692415   0.3399810435849;...
        0.0   0.0               0.0               0.8611363115941];
    WGT=[2.0, 1.0, 0.5555555555555, 0.3478548451375;...
        0.0, 1.0, 0.8888888888888, 0.6521451548625;...
        0.0, 0.0, 0.5555555555555, 0.6521451548625;...
        0.0, 0.0, 0.0,             0.3478548451375];
    
    D=E/(1+NU)/(1-2*NU)*[1-NU NU 0;NU 1-NU 0;0 0 (1-2*NU)/2];
    
    
    for LX=1:NINT  
        RI=XG(LX,NINT);
        for LY=1:NINT
            SI=XG(LY,NINT);
            
            [B,DET]=STDM(RI,SI,XX,B);
            
            WT = WGT(LX,NINT)*WGT(LY,NINT)*DET;
            S=S+B'*D*B*WT;
            
        end
    end
    for J=1:8      
        for I=J:8    
            S(J,I)=S(I,J);
        end
    end
%   SRC/Mechanics/ADDBAN.m
    ADDBAN(S, LM(:, N));
end    

% The third time stamp
cdata.TIM(3, :) = clock;

end




