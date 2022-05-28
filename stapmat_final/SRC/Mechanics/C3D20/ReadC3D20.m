% % - Called by
% %          C3D20Stiff.m
% % 
% % - Call procedures:
% %        ReadMaterial()
% %        ReadElements()

function ReadC3D20(NUMEG_ID)

% Read Material information
ReadMaterial()

% Read Element information
ReadElements(NUMEG_ID)

% the second time stamp
global cdata;
cdata.TIM(2,:) = clock;

end % end of ReadC3D20


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
fprintf(IOUT, '     SET          YOUNG''S       POSSION RATIO         DENSITY\n');
fprintf(IOUT, '   NUMBER         MODULUS\n');
fprintf(IOUT, '                     E                  nu              rho\n');

% Read material datas
sdata.NUME = cdata.NPAR(2);% 该单元组内的单元数
sdata.NUMMAT = cdata.NPAR(3);% 材料/截面属性数
NUMMAT = cdata.NPAR(3);
sdata.NGaussian = cdata.NPAR(4); % 该单元类型的高斯积分阶数
sdata.E = zeros(NUMMAT, 1, 'double');% 杨氏模量
sdata.nu = zeros(NUMMAT, 1, 'double');% 泊松比
sdata.rho = zeros(NUMMAT, 1, 'double');% 密度
for I = 1:cdata.NPAR(3)% 材料/截面属性数 
    tmp = str2num(fgetl(IIN));
    N = round(tmp(1));
    sdata.E(N) = tmp(2);
    sdata.nu(N) = tmp(3);
    sdata.rho(N) = tmp(4);
    fprintf(IOUT, ' %5d         %12.5e    %14.6e     %5d\n', N, tmp(2), tmp(3), tmp(4));
end

end % end of  ReadMaterial


% Read elements information
function ReadElements(NUMEG_ID)

global cdata;
global sdata;

% Get file pointer
IIN = cdata.IIN;
IOUT = cdata.IOUT;

fprintf(IOUT, '\n\n E L E M E N T   I N F O R M A T I O N\n');
fprintf(IOUT, '\n      ELEMENT \t ');
for i_out = 1:20
    fprintf(IOUT, 'NODE \t ');
end
fprintf(IOUT, '   MATERIAL\n');
fprintf(IOUT, '      NUMBER-N \t   ');
for i_out = 1:20
    fprintf(IOUT, 'I%d 	   ',i_out);
end
fprintf(IOUT, ' SET NUMBER\n');

% Get Position data
NUME = cdata.NPAR(2);
sdata.XYZ = zeros(sdata.NNODE * 3, NUME, 'double'); % 每一列为一个单元，每一列中的数表示该单元上的节点的xyz坐标
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

    if (NUMEG_ID ~= 1)
        for i3 = 1:sdata.NNODE
            sdata.ELEII(sdata.NUMEGEM(NUMEG_ID-1)+N,i3) = II(i3);
        end
    else
        for i3 = 1:sdata.NNODE
            sdata.ELEII(N,i3) = II(i3);
        end
    end
    for i3 = 1:sdata.NNODE
        sdata.ELEIIC3D20(N,i3) = II(i3);
    end

    MTYPE = round(tmp(sdata.NNODE+2));

    %   Save element information
    for i1 = 1:sdata.NNODE
        XYZ( i1*3 - 2, N ) = X( II( i1 ) );
        XYZ( i1*3 - 1, N ) = Y( II( i1 ) );
        XYZ( i1*3, N ) = Z( II( i1 ) );
    end
    
    MATP(N) = MTYPE;

    fprintf(IOUT, '     %8d',N);
    for i_out = 1:20
        fprintf(IOUT, '%8d',II(i_out));
    end
    fprintf(IOUT, '     %5d\n',MTYPE);
    
    % Compute connectivity matrix
    for i2 = 1:sdata.NNODE
        LM( i2*3-2, N ) = ID( 1, II(i2) );
        LM( i2*3-1, N ) = ID( 2, II(i2) );
        LM( i2*3  , N ) = ID( 3, II(i2) );
    end

    %   Updata column heights and bandwidth
    ColHt(LM(:, N))
end

sdata.XYZ = XYZ; sdata.MATP = MATP; sdata.LM = LM;


end % end of ReadElements