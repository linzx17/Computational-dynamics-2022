% ******************************************************************
%*                                                                 *
%*                         S T A P M A T                           *
%*                                                                 *
%*      An In-CORE SOLUTION STATIC ANALYSIS PROGRAM IN MATLAB      *
%*      Adapted from STAP90 (FORTRAN 90) for teaching purpose      *
%*                                                                 *
%*  Computational Dynamics Group, School of Aerospace Engineering  *
%*  Tsinghua University, 2019.02.20                                *
%*                                                                 *
%* *****************************************************************

% Set paths of functions
AddPath();

% Define Global Variables
global cdata;
global sdata;
cdata = ControlData;
sdata = SolutionData;

% Read InPut file
fileID = 'Job-C3D8-cantilever-beam-rho7800-dyna';
% fileID = 'C3D20_cantilever_beam';
% fileID = 'Job_C3D20-cantilever-beam-element320';
% fileID = 'C3D8_3surf';
% mkdir Data\ C3D20_cantilever_beam;

fname = [fileID,'.in'];              % Specify the file name

ReadFile(fname); % 读取标题行、控制行、节点数据、载荷数

% Write basic data of program 

cdata.IOUT = fopen(['.\Data\',fileID,'\',fileID,'.OUT'], 'w');
WriteParasOut(); % 创建输出文件

cdata.IDAT_ANIM = fopen(['.\Data\',fileID,'\',fileID,'_anim.DAT'], 'w');
WriteParasDatANIM();%创建Tecplot可读文件

cdata.IDAT_CURV = fopen(['.\Data\',fileID,'\',fileID,'_curv.DAT'], 'w');
WriteParasDatCURV();%创建Tecplot可读文件

% Form the stiffness matrix
GetStiff();

% % solve
SOLVE();

% % Triangularize stiffness matrix
% Solve();
% 
% % 特征值和特征向量
% Eigenvalue();

% Finalize
Finalize();

% ----------------------- Functions -----------------------------------

% Functions
% Add paths of functions
function AddPath()
close all;
clear;
clc;

addpath .\SRC\Initiation
addpath .\SRC\BasicData
addpath .\SRC\Mechanics
addpath .\SRC\Mechanics\Truss
addpath .\SRC\Mechanics\Plane
addpath .\SRC\Mechanics\C3D8
addpath .\SRC\Mechanics\C3D20
addpath .\SRC\Solver
addpath .\SRC\Solver\Statics
addpath .\SRC\\Solver\Dynamics

end

function Finalize()
global cdata;
global sdata;
TIM = cdata.TIM;
time = zeros(5, 1, 'double');
time(1) = etime(TIM(2,:), TIM(1,:));
time(2) = etime(TIM(3,:), TIM(2,:));
time(3) = etime(TIM(4,:), TIM(3,:));
time(4) = etime(TIM(5,:), TIM(4,:));
time(5) = etime(TIM(6,:), TIM(5,:));
time(6) = etime(TIM(6,:), TIM(1,:));

fprintf(cdata.IOUT, ['\n\n' ...
    ' S O L U T I O N   T I M E   L O G   I N   S E C\n\n' ...
    '     TIME FOR INPUT PHASE  . . . . . . . . . . . . . . = %12.2f\n' ...
    '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . = %12.2f\n' ...
    '     TIME FOR FACTORIZATION OF STIFFNESS MATRIX  . . . = %12.2f\n' ...
    '     TIME FOR LOAD CASE SOLUTIONS  . . . . . . . . . . = %12.2f\n\n' ...
    '     TIME FOR Eigenvalues  . . . . . . . . . . . . . . = %12.2f\n\n' ...
    '      T O T A L   S O L U T I O N   T I M E  . . . . . = %12.2f\n'], ...
    time(1), time(2), time(3), time(4), time(5), time(6));

fprintf(['\n' ...
    ' S O L U T I O N   T I M E   L O G   I N   S E C\n\n' ...
    '     TIME FOR INPUT PHASE  . . . . . . . . . . . . . . = %12.2f\n' ...
    '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . = %12.2f\n' ...
    '     TIME FOR FACTORIZATION OF STIFFNESS MATRIX  . . . = %12.2f\n' ...
    '     TIME FOR LOAD CASE SOLUTIONS  . . . . . . . . . . = %12.2f\n' ...
    '     TIME FOR Eigenvalues  . . . . . . . . . . . . . . = %12.2f\n\n' ...
    '      T O T A L   S O L U T I O N   T I M E  . . . . . = %12.2f\n'], ...
    time(1), time(2), time(3), time(4), time(5), time(6));

for N = 1:cdata.NUMEM
    for NN = 1:sdata.NNODE
        fprintf(cdata.IDAT_ANIM,'%7d',sdata.ELEII(cdata.NUMEM,NN));
    end
    fprintf(cdata.IDAT_ANIM,'\n');
end

fclose(cdata.IIN);
fclose(cdata.IOUT);
fclose(cdata.IDAT_ANIM);
fclose(cdata.IDAT_CURV);
end