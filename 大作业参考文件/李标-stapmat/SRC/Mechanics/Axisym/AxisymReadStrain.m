%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read the element information of truss                       *
%*                                                                 *
%* - Call procedures:                                              *
%*     ReadTruss.m - ReadMaterial()                                *
%*     ReadTruss.m - ReadElements()                                *
%*                                                                 *
%* - Called by :                                                   *
%*     ./StrainStiff.m                                             *
%*                                                                 *
%* *****************************************************************
%增添了网格输出
function AxisymReadStrain()

% Read Material information
ReadMaterial()

% Read Element information
ReadElements()

% the second time stamp
global cdata;
cdata.TIM(2,:) = clock;

end

% ----------------------- Functions -----------------------------------
% Read Material information
function ReadMaterial()

global cdata;
global sdata;
% Get file pointers
IIN = cdata.IIN;
IOUT = cdata.IOUT;

if (cdata.NPAR(3) == 0) 
    cdata.NPAR(3) = 1; 
end
fprintf(IOUT, '\n M A T E R I A L   D E F I N I T I O N\n');
fprintf(IOUT, '\n NUMBER OF DIFFERENT SETS OF MATERIAL\n');
fprintf(IOUT, ' AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . = %10d\n', ...
    cdata.NPAR(3));
fprintf(IOUT, '  SET       YOUNG''S       POISSON''S         GAUSS NUMERICAL  \n');
fprintf(IOUT, ' NUMBER     MODULUS         RATIO       INTEGRATION ORDER \n');
fprintf(IOUT, '               E              PR                  NINT \n');


% Read material datas
sdata.NUME = cdata.NPAR(2);
sdata.NUMMAT = cdata.NPAR(3);
NUMMAT = cdata.NPAR(3);
sdata.E = zeros(NUMMAT, 1, 'double');
sdata.AREA = zeros(NUMMAT, 1, 'double');

for I = 1:cdata.NPAR(3)
    tmp = str2num(fgetl(IIN));
    N = round(tmp(1));  %最近的整数
    sdata.E(N) = tmp(2);
    sdata.PR(N) = tmp(3);
    sdata.density(N) = tmp(4);
    sdata.NINT(N) = tmp(5);  %将平面单元的厚度删除，仅仅保存高斯积分点，这里用的是减缩积分
    sdata.masschoose = round(tmp(6)); %质量矩阵的生成方式。
    sdata.radius_way = round(tmp(7)); %质量矩阵的生成方式。
 %下面的fprintf也就进行了修改，仅仅输出上面的几个元素，将tmp（5）删除了
    fprintf(IOUT, '%5d    %12.5e        %2.1f       %14.6e \n', N, tmp(2), tmp(3), tmp(4));
end

end   %

% Read elements information
function ReadElements()

global cdata;
global sdata;

% Get file pointer
IIN = cdata.IIN;
IOUT = cdata.IOUT;

fprintf(IOUT, '\n\n E L E M E N T   I N F O R M A T I O N\n'); %在原来的基础上添加了4个节点
fprintf(IOUT, '\n      ELEMENT          NODE        NODE        NODE        NODE              NODE        NODE        NODE        NODE      MATERIAL\n');
fprintf(IOUT, '      NUMBER-N           I           J           K           L          A          B           C         D       SET NUMBER\n');

% Get Position data
NUME = cdata.NPAR(2);
sdata.XYZ = zeros(16, NUME, 'double');   %将8改成了16
sdata.MATP = zeros(NUME, 1, 'int64');                 % the type of material
sdata.LM = zeros(16, NUME, 'int64');                 %将8改成了16 % connectivity matrix
sdata.MHT = zeros(sdata.NEQ, 1, 'int64');
X = sdata.X; Y = sdata.Y; Z = sdata.Z; ID = sdata.ID;
XYZ = sdata.XYZ; MATP = sdata.MATP; LM = sdata.LM;
NUMNP = cdata.NUMNP;
filename = strcat('.\Data\',cdata.HED);
mkdir(filename);
filename = strcat(filename,'\mid.vtk');
fid = fopen(filename, 'w'); 
fprintf(fid,[ '# vtk DataFile Version 3.0\n'...
              'VTK from Matlab\n'...
              'ASCII\n'...
              'DATASET POLYDATA\n']);
fprintf(fid,'POINTS	%d	double\n',NUMNP);
for i = 1:NUMNP
    fprintf(fid,'%e %e %e\n',X(i),Y(i),Z(i));
end
fprintf(fid,'POLYGONS %d %d\n',NUME,9*NUME);

for N = 1:NUME
    tmp = str2num(fgetl(IIN));
    I = round(tmp(2));
    J = round(tmp(3));
    K = round(tmp(4));
    L = round(tmp(5)); 
    A = round(tmp(6));%增添了4个节点,ABCD
    B = round(tmp(7));
    C = round(tmp(8));
    D = round(tmp(9));
    MTYPE = round(tmp(10)); %将6改成了10
%   Save element information
    XYZ(1, N) = X(I);
    XYZ(2, N) = Y(I);
    XYZ(3, N) = X(J);
    XYZ(4, N) = Y(J);
    XYZ(5, N) = X(K);
    XYZ(6, N) = Y(K);
    XYZ(7, N) = X(L);
    XYZ(8, N) = Y(L);
    XYZ(9, N) = X(A);   %从原来的8添加成了16
    XYZ(10, N) = Y(A);
    XYZ(11, N) = X(B);
    XYZ(12, N) = Y(B);
    XYZ(13, N) = X(C);
    XYZ(14, N) = Y(C);
    XYZ(15, N) = X(D);
    XYZ(16, N) = Y(D);
    MATP(N) = MTYPE;
    %下面的fprintf也就行了修改
    fprintf(IOUT, '%10d      %10d  %10d  %10d  %10d %10d  %10d  %10d  %10d       %5d\n', N, I, J, K, L, A, B, C, D  ,MTYPE);
    fprintf(fid,'8 %d %d %d %d %d %d %d %d\n',I-1, A-1, J-1, B-1, K-1, C-1, L-1, D-1);
%   Compute connectivity matrix
    LM(1, N) = ID(1, I);
    LM(3, N) = ID(1, J);
    LM(5, N) = ID(1, K);
    LM(7, N) = ID(1, L);
    LM(9, N) = ID(1, A);   %ID的第一行是X方向的自由度对应的整体刚度矩阵的位置，增了4个
    LM(11, N) = ID(1, B);
    LM(13, N) = ID(1, C);
    LM(15, N) = ID(1, D);
    
    LM(2, N) = ID(2, I);
    LM(4, N) = ID(2, J);
    LM(6, N) = ID(2, K);
    LM(8, N) = ID(2, L);  %ID的第二行是Y方向的自由度对应的整体刚度矩阵的位置，增了4个
    LM(10, N) = ID(2, A);
    LM(12, N) = ID(2, B);
    LM(14, N) = ID(2, C);
    LM(16, N) = ID(2, D); 

end
sdata.XYZ = XYZ; sdata.MATP = MATP; sdata.LM = LM;
sdata.X = double(0);
sdata.Y = double(0);
sdata.Z = double(0);

end