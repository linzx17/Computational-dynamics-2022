% % - Called by
% %          PlaneStiff.m
% % 
% % - Call procedures:
% %        ReadMaterial()
% %        ReadElements()


function ReadPlane()

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
fprintf(IOUT, '     SET          YOUNG''S            POSSION         NUMBER OF GAUSSIAN\n');
fprintf(IOUT, ' NUMBER     MODULUS             RATIO             INTEGRATION POINT   \n');
fprintf(IOUT, '                           E                       nu\n');


% Read material datas
sdata.NUME = cdata.NPAR(2);% 该单元组内的单元数
sdata.NUMMAT = cdata.NPAR(3);% 材料/截面属性数
NUMMAT = cdata.NPAR(3);
sdata.E = zeros(NUMMAT, 1, 'double');% 杨氏模量
sdata.nu = zeros(NUMMAT, 1, 'double');% 泊松比
sdata.NGaussian = zeros(NUMMAT, 1, 'int64');% 高斯积分点个数
for I = 1:cdata.NPAR(3)% 材料/截面属性数 
    tmp = str2num(fgetl(IIN));
    N = round(tmp(1));
    sdata.E(N) = tmp(2);
    sdata.nu(N) = tmp(3);
    sdata.NGaussian(N) = tmp(4);
    fprintf(IOUT, '  %5d          %12.5e    %14.6e                   %5d\n', N, tmp(2), tmp(3), tmp(4));
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
fprintf(IOUT, '\n      ELEMENT \t NODE \t NODE \t NODE \t NODE \t   MATERIAL\n');
fprintf(IOUT, '      NUMBER-N \t   I1 \t    I2 \t    I3 \t    I4 \t SET NUMBER\n');

% Get Position data
NUME = cdata.NPAR(2);
sdata.XYZ = zeros(sdata.NNODE * 3, NUME, 'double'); % 每一列为一个单元，12行为4个节点上12个xyz坐标
sdata.MATP = zeros(NUME, 1, 'int64');                 % the type of material
sdata.LM = zeros(sdata.NNODE * sdata.NDOF, NUME, 'double');                  % connectivity matrix
sdata.MHT = zeros(sdata.NEQ, 1, 'int64');
X = sdata.X; Y = sdata.Y; Z = sdata.Z; ID = sdata.ID;
XYZ = sdata.XYZ; MATP = sdata.MATP; LM = sdata.LM;

for N = 1:NUME
    tmp = str2num(fgetl(IIN));

    II = zeros(sdata.NNODE,1); % 该单元上的节点编号
    for ii = 1:sdata.NNODE
        II(ii) = round(tmp(ii+1));
    end

    MTYPE = round(tmp(sdata.NNODE+2));
%     I1 = round(tmp(2));
%     I2 = round(tmp(3));
%     I3 = round(tmp(4));
%     I4 = round(tmp(5));
%     MTYPE = round(tmp(6));
    
%   Save element information
    for i1 = 1:sdata.NNODE
        XYZ( i1*3 - 2, N ) = X( II( i1 ) );
        XYZ( i1*3 - 1, N ) = Y( II( i1 ) );
        XYZ( i1*3, N ) = Z( II( i1 ) );
    end
%     XYZ(1, N) = X(I1);
%     XYZ(2, N) = Y(I1);
%     XYZ(3, N) = Z(I1);
%     XYZ(4, N) = X(I2);
%     XYZ(5, N) = Y(I2);
%     XYZ(6, N) = Z(I2);
%     XYZ(7, N) = X(I3);
%     XYZ(8, N) = Y(I3);
%     XYZ(9, N) = Z(I3);
%     XYZ(10, N) = X(I4);
%     XYZ(11, N) = Y(I4);
%     XYZ(12, N) = Z(I4);

    MATP(N) = MTYPE;
    
    fprintf(IOUT, '   %10d            %10d      %10d      %10d      %10d           %5d\n', N, II(1), II(2), II(3), II(4), MTYPE);

%   Compute connectivity matrix
    for i2 = 1:sdata.NNODE
        LM( i2*2-1, N ) = ID( 1, II(i2) );
        LM( i2*2, N) = ID( 2, II(i2) );
    end
%     LM(1, N) = ID(1, I1);
%     LM(2, N) = ID(2, I1);
%     LM(3, N) = ID(1, I2);
%     LM(4, N) = ID(2, I2);
%     LM(5, N) = ID(1, I3);
%     LM(6, N) = ID(2, I3);
%     LM(7, N) = ID(1, I4);
%     LM(8, N) = ID(2, I4);

%   Updata column heights and bandwidth
    ColHt(LM(:, N))
end
sdata.XYZ = XYZ; sdata.MATP = MATP; sdata.LM = LM;

% Clear the memory of X, Y, Z
sdata.X = double(0);
sdata.Y = double(0);
sdata.Z = double(0);

end