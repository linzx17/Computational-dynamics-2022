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

function planeReadStrain()

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
fprintf(IOUT, '  SET       YOUNG''S       POISSON''S         ELEMENT     GAUSS NUMERICAL  \n');
fprintf(IOUT, ' NUMBER     MODULUS         RATIO          THICKNESS   INTEGRATION ORDER \n');
fprintf(IOUT, '               E              PR              THIC            NINT \n');


% Read material datas
sdata.NUME = cdata.NPAR(2);
sdata.NUMMAT = cdata.NPAR(3);
NUMMAT = cdata.NPAR(3);
sdata.E = zeros(NUMMAT, 1, 'double');
sdata.AREA = zeros(NUMMAT, 1, 'double');
for I = 1:cdata.NPAR(3)
    tmp = str2num(fgetl(IIN));
    N = round(tmp(1));
    sdata.E(N) = tmp(2);
    sdata.PR(N) = tmp(3);
    sdata.THIC(N) = tmp(4);
    sdata.NINT(N) = tmp(5);
    fprintf(IOUT, '%5d    %12.5e        %2.1f       %14.6e        %2d\n', N, tmp(2), tmp(3), tmp(4), tmp(5));
end

end

% Read elements information
function ReadElements()

global cdata;
global sdata;

% Get file pointer
IIN = cdata.IIN;
IOUT = cdata.IOUT;

fprintf(IOUT, '\n\n E L E M E N T   I N F O R M A T I O N\n');
fprintf(IOUT, '\n      ELEMENT          NODE        NODE        NODE        NODE       MATERIAL\n');
fprintf(IOUT, '      NUMBER-N           I           J           K           L       SET NUMBER\n');

% Get Position data
NUME = cdata.NPAR(2);
sdata.XYZ = zeros(8, NUME, 'double');
sdata.MATP = zeros(NUME, 1, 'int64');                 % the type of material
sdata.LM = zeros(8, NUME, 'int64');                  % connectivity matrix
sdata.MHT = zeros(sdata.NEQ, 1, 'int64');
X = sdata.X; Y = sdata.Y; ID = sdata.ID;
XYZ = sdata.XYZ; MATP = sdata.MATP; LM = sdata.LM;


for N = 1:NUME
    tmp = str2num(fgetl(IIN));
    I = round(tmp(2));
    J = round(tmp(3));
    L = round(tmp(4));
    K = round(tmp(5));
    MTYPE = round(tmp(6));
    
%   Save element information
    XYZ(1, N) = X(I);
    XYZ(2, N) = Y(I);
    XYZ(3, N) = X(J);
    XYZ(4, N) = Y(J);
    XYZ(5, N) = X(K);
    XYZ(6, N) = Y(K);
    XYZ(7, N) = X(L);
    XYZ(8, N) = Y(L);
    MATP(N) = MTYPE;
    
    fprintf(IOUT, '%10d      %10d  %10d  %10d  %10d       %5d\n', N, I, J, K, L, MTYPE);

%   Compute connectivity matrix
    LM(1, N) = ID(1, I);
    LM(3, N) = ID(1, J);
    LM(5, N) = ID(1, K);
    LM(7, N) = ID(1, L);
    LM(2, N) = ID(2, I);
    LM(4, N) = ID(2, J);
    LM(6, N) = ID(2, K);
    LM(8, N) = ID(2, L);
%   Updata column heights and bandwidth
    ColHt(LM(:, N))
end
sdata.XYZ = XYZ; sdata.MATP = MATP; sdata.LM = LM;

% Clear the memory of X, Y, Z
sdata.X = double(0);
sdata.Y = double(0);
sdata.Z = double(0);

end