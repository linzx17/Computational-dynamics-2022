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
fname = strcat('.\Data\', fname);           % Deal the filename

% Get global class
global cdata;
global sdata;

% Open files
cdata.IIN = fopen(fname, 'r');

% Begin Read input file
fprintf('Input phase ...\n\n');

% the first time stamp
cdata.TIM = zeros(4, 6, 'double');
cdata.TIM(1,:) = clock;

IIN = cdata.IIN;
%% Read Control data
cdata.HED = fgetl(IIN);

tmp = str2num(fgetl(IIN));
cdata.NUMNP = int64(tmp(1));
cdata.NUMEG = int64(tmp(2));
cdata.NLCASE = int64(tmp(3));
cdata.MODEX = int64(tmp(4));

if (cdata.NUMNP == 0) 
    return; 
end


    tmp = str2num(fgetl(IIN));
    cdata.dt = double(tmp(1));
    cdata.nt = double(tmp(2));
    cdata.gamma = double(tmp(3));
    cdata.beta = double(tmp(4));
    cdata.betamode = int64(tmp(5));
    cdata.outt = int64(tmp(6));
    if (cdata.betamode == 1 && cdata.gamma ~= 0.5)
        error(' *** ERROR *** gamma~=0.5,不可参数优化');
    end


%% Read nodal point data
InitBasicData();
% Define local variables to speed
ID = sdata.ID; X = sdata.X; Y = sdata.Y; Z = sdata.Z;
for i = 1:cdata.NUMNP
    tmp = str2num(fgetl(IIN));
    ID(1, tmp(1)) = int64(tmp(2));
    ID(2, tmp(1)) = int64(tmp(3));
    ID(3, tmp(1)) = int64(tmp(4));
    X(tmp(1)) = double(tmp(5));
    Y(tmp(1)) = double(tmp(6));
    Z(tmp(1)) = double(tmp(7));
end
sdata.ID = ID; sdata.X = X; sdata.Y = Y; sdata.Z = Z;
%% Compute the number of equations
sdata.IDOrigin = ID;
NEQ = 0;
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
%% Read load data
% Init control data
NLCASE = cdata.NLCASE;
sdata.R = zeros(NEQ, NLCASE, 'double');
R = sdata.R;
% Read data
for N = 1:cdata.NLCASE
    tmp = str2num(fgetl(IIN));
    cdata.LL = int64(tmp(1)); cdata.NLOAD = int64(tmp(2));
    NLOAD = cdata.NLOAD;
%   Init load data
    sdata.NOD = zeros(NLOAD, 1, 'int64');
    sdata.IDIRN = zeros(NLOAD, 1, 'int64');
    sdata.FLOAD = zeros(NLOAD, 1, 'double');
    sdata.CHNOD = zeros(NLOAD, 3, 'double');
    NOD = sdata.NOD; IDIRN = sdata.IDIRN; FLOAD = sdata.FLOAD;
    CHNOD = sdata.CHNOD;
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
    CHNOD(:,2) = X(NOD);
    CHNOD(:,3) = Y(NOD);
%   Compute load vector
    for L = 1:NLOAD
        II = ID(IDIRN(L), NOD(L));
        CHNOD(L,1) = II;
        if (II > 0) 
            R(II, N) = R(II, N) + FLOAD(L); 
        end
    end
    sdata.NOD = NOD; sdata.IDIRN = IDIRN; sdata.FLOAD = FLOAD; sdata.R = R;
    sdata.CHNOD = CHNOD;
end

end

%% Functions


% InitBasicData
function InitBasicData()
global cdata;
global sdata;

cdata.NPAR = zeros(10, 1, 'int64');

sdata.ID = zeros(3,cdata.NUMNP, 'int64');
sdata.X = zeros(cdata.NUMNP, 1, 'double');
sdata.Y = zeros(cdata.NUMNP, 1, 'double');
sdata.Z = zeros(cdata.NUMNP, 1, 'double');
end