%* *****************************************************************
%* - Function of STAPMAT in initialization phase                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Read input file of STAPMAT                                  *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Initiation/ReadFile.m - InitBasicData()                 *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************

function ReadFile(fname)
fname = strcat('.\Data\', fname);           % Deal the filename % strcat() 水平串联字符串

% Get global class
global cdata;
global sdata;

% Open files
cdata.IIN = fopen(fname, 'r'); % 文件标识符（>=3的整数）

% Begin Read input file
fprintf('Input phase ...\n\n');

% the first time stamp
cdata.TIM = zeros(6, 6, 'double');
cdata.TIM(1,:) = clock;

IIN = cdata.IIN;
%% Read Control data
cdata.HED = fgetl(IIN); % 读取文件下一行并删除换行符

tmp = str2num(fgetl(IIN));
cdata.NUMNP = int64(tmp(1)); % 节点总数
cdata.NUMEG = int64(tmp(2)); % 单元组总数
cdata.NUMEM = int64(tmp(3)); % 单元总数
cdata.NLCASE = int64(tmp(4)); % 载荷工况数
cdata.MODEX = int64(tmp(5)); % 求解模式

if cdata.MODEX == 4 % 如果求解模式为4，即求解动力学问题，则读求解时间和求解方法
    cdata.DSTIME = tmp(6);
    cdata.DSMETH = int64(tmp(7));
end

if (cdata.NUMNP == 0)
    return;
end

%% Read nodal point data
InitBasicData();
% Define local variables to speed
ID = sdata.ID; X = sdata.X; Y = sdata.Y; Z = sdata.Z;
for i = 1:cdata.NUMNP
    tmp = str2num(fgetl(IIN));
    ID(1, i) = int64(tmp(2));
    ID(2, i) = int64(tmp(3));
    ID(3, i) = int64(tmp(4));
    X(i) = double(tmp(5));
    Y(i) = double(tmp(6));
    Z(i) = double(tmp(7));
end
sdata.ID = ID; sdata.X = X; sdata.Y = Y; sdata.Z = Z;
%% Compute the number of equations
sdata.IDOrigin = ID;
NEQ = 0; % 自由度个数
for N=1:cdata.NUMNP
    for I=1:3
        if (ID(I,N) == 0)
            NEQ = NEQ + 1;
            ID(I,N) = NEQ;
        else
            ID(I,N) = 0;
        end
    end
end
sdata.ID = ID;
sdata.NEQ = NEQ;
sdata.MASS_DIRCTION = zeros(NEQ,NEQ);
%% Read load data
% Init control data
NLCASE = cdata.NLCASE; % 载荷工况数
sdata.R = zeros(NEQ, NLCASE, 'double'); % R(i,j)为第i个自由度上第j个载荷工况的载荷值
R = sdata.R;
% Read data
for N = 1:cdata.NLCASE % 载荷工况数
    tmp = str2num(fgetl(IIN));
    cdata.LL = int64(tmp(1)); % 载荷工况号
    cdata.NLOAD = int64(tmp(2)); % 本工况中集中载荷的个数
% % % % 如果求解动力学问题，则存储第一个载荷工况下的动载荷类型、动载荷参数
    if (4 == cdata.MODEX) && (1 == N)
        cdata.DLTYPE = int64(tmp(3));
        cdata.DLPAR = tmp(4:7);
    end
% % % % 
    NLOAD = cdata.NLOAD;
%   Init load data
    sdata.NOD = zeros(NLOAD, 1, 'int64'); % 集中载荷作用的节点号
    sdata.IDIRN = zeros(NLOAD, 1, 'int64'); % 载荷作用方向
    sdata.FLOAD = zeros(NLOAD, 1, 'double'); % 载荷作用方向（1-x方向，2-y方向，3-z方向）
    NOD = sdata.NOD; IDIRN = sdata.IDIRN; FLOAD = sdata.FLOAD;
    
%   Read load data
    for I = 1:NLOAD
        tmp = str2num(fgetl(IIN));
        NOD(I) = int64(tmp(1));
        IDIRN(I) = int64(tmp(2));
        FLOAD(I) = double(tmp(3));
    end
    if (cdata.MODEX == 0) 
        return;
    end
    
%   Compute load vector
    for L = 1:NLOAD
        II = ID(IDIRN(L), NOD(L));
        if (II > 0)
            R(II, N) = R(II, N) + FLOAD(L);
        end
    end
    sdata.NOD = NOD;
    sdata.IDIRN = IDIRN;
    sdata.FLOAD = FLOAD;
    sdata.R = R;
end

end

%% Functions


% InitBasicData
function InitBasicData()
global cdata;
global sdata;

cdata.NPAR = zeros(10, 1, 'int64');

sdata.ID = zeros(3,cdata.NUMNP, 'int64'); 
% ID数组的第1、2、3行分别表示xyz方向自由度，第j(1<=j<=cdata.NUMNP)列表示第j个节点
% ID(i,j)的值表示第j个节点的第i个自由度是否使用(0=free, 1=deleted)
sdata.X = zeros(cdata.NUMNP, 1, 'double'); % 节点x坐标
sdata.Y = zeros(cdata.NUMNP, 1, 'double'); % 节点y坐标
sdata.Z = zeros(cdata.NUMNP, 1, 'double'); % 节点z坐标
end