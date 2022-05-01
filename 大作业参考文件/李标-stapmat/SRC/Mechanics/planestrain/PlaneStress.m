%* *****************************************************************
%* - Function of STAPMAT in solver phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     To calculate stresses                                       *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Solver/GetStress.m                                      *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.22                *
%*                                                                 *
%* *****************************************************************

function PlaneStress(NUM,NG)

% Get global data
global cdata;
global sdata;
IOUT = cdata.IOUT;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
EN = sdata.E; LM = sdata.LM; NNU = sdata.PR;
U = sdata.DIS(:, NUM);
XX = zeros(2, 4, 'double');
UV = zeros(1, 8, 'double');
B = zeros(3, 8, 'double');
maxx = 0.0;
maxxn = 0;
fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '       ELEMENT             SoildSX            SoildSY'...
    '            SoildSXY\n' ...
    '       NUMBER\n'], NG);


for N = 1:NUME
    MTYPE = MATP(N);
    NU = NNU(MTYPE);
    E = EN(MTYPE);
    D=E/(1+NU)/(1-2*NU)*[1-NU NU 0;NU 1-NU 0;0 0 (1-2*NU)/2];
    UV (1,1:8) =0;
    for i = 1:2
        for j = 1:4
            XX(i,j) = XYZ(2*j+i-2,N);
        end
    end
    
    for i = 1:8
        if LM(i, N) > 0
        UV(1,i) = U(LM(i, N));
        end
    end
    [B1,~]=STDM(1,1,XX,B);
    [B2,~]=STDM(-1,1,XX,B);
    [B3,~]=STDM(1,-1,XX,B);
    [B4,~]=STDM(-1,-1,XX,B);
%     [B1,~]=STDM(3^0.5/3,3^0.5/3,XX,B);
%     [B2,~]=STDM(-3^0.5/3,3^0.5/3,XX,B);
%     [B3,~]=STDM(3^0.5/3,-3^0.5/3,XX,B);
%     [B4,~]=STDM(-3^0.5/3,-3^0.5/3,XX,B);
    STR1 = D*(B1*UV');
    STR2 = D*(B2*UV');
    STR3 = D*(B3*UV');
    STR4 = D*(B4*UV');
    STR = (STR1+STR2+STR3+STR4)/4;
    fprintf(IOUT, ' %10d           %13.6e    %13.6e    %13.6e\n'...
        , N, STR(1), STR(2), STR(3));
    POINTSX(1) = STR1(1);
    POINTSX(2) = STR2(1);
    POINTSX(3) = STR3(1);
    POINTSX(4) = STR4(1);
    for ii = 1:4
        if(POINTSX(ii)>maxx)%最大x方向应力用以验证小孔应力集中
            maxx = max(POINTSX);
            maxcenx=XYZ(2*ii-1,N);
            maxceny=XYZ(2*ii,N);
        end
    end
end

fprintf(IOUT, ['\n\n  MAX SolidSX  \n\n' ... 
'       SolidSX            节点坐标\n\n']);
fprintf(IOUT, '    %13.6e      (%5.4f , %5.4f) \n',...
     maxx, maxcenx, maxceny);
    