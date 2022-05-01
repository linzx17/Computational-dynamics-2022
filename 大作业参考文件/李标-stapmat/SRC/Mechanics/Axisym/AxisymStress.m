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

function AxisymStress(NUM,NG)

% Get global data
global cdata;
global sdata;
IOUT = cdata.IOUT;
MODEX = cdata.MODEX;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
EN = sdata.E; LM = sdata.LM; NNU = sdata.PR;
U = sdata.DIS(:, NUM);
sdata.STRAIN = zeros(NUME, 4, 'double');
sdata.STRESS = zeros(NUME, 5, 'double');
STRAIN = sdata.STRAIN;
STRESS = sdata.STRESS;
radius_way = sdata.radius_way;
XX = zeros(2, 8, 'double');   %将4改成了8
UV = zeros(1, 16, 'double');  %这应该是用来储存位移的一个中介量，将8改成了16
 %轴对称单元并且是轴对称载荷
%输出增添了轴对称的应力
if (MODEX == 1)
    fprintf(IOUT, ['\n\n  S T R E S S  C A L C U L A T I O N S  F O R  ' ...
        'E L E M E N T  G R O U P %4d\n\n' ...
        '       ELEMENT             SoildSR            SoildSZ'...
        '            SoildSRZ            SoildStheta  \n' ...
        '       NUMBER\n'], NG);
end
MTYPE = MATP(1);
NU = NNU(MTYPE);
E = EN(MTYPE);
D=E/(1+NU)/(1-2*NU)*[1-NU NU 0 NU;NU 1-NU 0 NU;0 0 (1-2*NU)/2 0;NU NU 0 1-NU]; %改成了轴对称单元的D
for N = 1:NUME
    UV (1,1:16) =0;%将8改成了16
    for i = 1:2
        for j = 1:8 %将4改成了8
            XX(i,j) = XYZ(2*j+i-2,N);
        end
    end
    for i = 1:16 %将8改成了16
        if LM(i, N) > 0
            UV(1,i) = U(LM(i, N));
        end
    end
    [B1, ~, ~] = generateBAD (1, 1, XX, radius_way);
    [B2, ~, ~] = generateBAD (2, 1, XX, radius_way);
    [B3, ~, ~] = generateBAD (2, 2, XX, radius_way);
    [B4, ~, ~] = generateBAD (1, 1, XX, radius_way);
    STI1 = (B1*UV');
    STI2 = (B2*UV');
    STI3 = (B3*UV');
    STI4 = (B4*UV');
    
    STR1 = D*STI1; %每个点有4个应力
    STR2 = D*STI2;
    STR3 = D*STI3;
    STR4 = D*STI4;
    STI = (STI1+STI2+STI3+STI4)/4;
    STR = (STR1+STR2+STR3+STR4)/4;
    Mises = (((STR(1)-STR(2))^2+(STR(2)-STR(4))^2+(STR(1)-STR(4))^2+6*STR(3)^2)/2)^0.5;
    if (MODEX == 1)
        fprintf(IOUT, ' %10d           %13.6e    %13.6e    %13.6e    %13.6e\n'...
            , N, STR(1), STR(2), STR(3), STR(4));
        %增添了高斯点的应力输出，输出从左下角逆时针循环开始
        fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
            , STR1(1), STR1(2), STR1(3), STR1(4));
        fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
            , STR2(1), STR2(2), STR2(3), STR2(4));
        fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
            , STR3(1), STR3(2), STR3(3), STR3(4));
        fprintf(IOUT, '                %13.6e    %13.6e    %13.6e    %13.6e\n'...
            , STR4(1), STR4(2), STR4(3), STR4(4));
    end
        STRAIN(N,:) = STI;
        STRESS(N,1:4) = STR;
        STRESS(N,5) = Mises;
end
[maxMises,index] = max(STRESS(:,5));
sdata.maxMises(1) = maxMises;
sdata.maxMises(2) = index;
sdata.STRAIN = STRAIN;
sdata.STRESS = STRESS;
end

