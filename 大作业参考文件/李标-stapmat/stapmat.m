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
% 直接运行为底部固定圆筒的外侧激励载荷，更改加载模式需手动调整参数
% Set paths of functions
function stapmat()
clear;
close all;
clc;
AddPath();

% Define Global Variables
global cdata;
global sdata;
cdata = ControlData;
sdata = SolutionData;

% Read InPut file
fname = 'A1xisym-unipressure(1)';              % Specify the file name
ReadFile(fname);
% Write basic data of program 
WriteParasOut(fname);
cdata.TIM(2, :) = clock;
% Form the stiffness matrix
 GetStiff();
cdata.TIM(3, :) = clock;
% Triangularize stiffness matrix
Solve();
cdata.TIM(4,:) = clock;
% Finalize
Finalize();

% ----------------------- Functions -----------------------------------

% Functions
% Add paths of functions
function AddPath()
clear;
close all;
clc;

addpath .\SRC\Initiation
addpath .\SRC\BasicData
addpath .\SRC\Mechanics
addpath .\SRC\Mechanics\Truss
addpath .\SRC\Mechanics\planestrain
addpath .\SRC\Mechanics\Axisym
addpath .\SRC\Solver
end

function Finalize()

TIM = cdata.TIM;
time = zeros(5, 1, 'double');
time(1) = etime(TIM(2,:), TIM(1,:));
time(2) = etime(TIM(3,:), TIM(2,:));
time(3) = etime(TIM(4,:), TIM(3,:));
time(4) = etime(TIM(4,:), TIM(1,:));

fprintf(cdata.IOUT, ['\n\n' ...
    ' S O L U T I O N   T I M E   L O G   I N   S E C\n\n' ...
    '     TIME FOR INPUT PHASE  . . . . . . . . . . . . . . = %12.2f\n' ...
    '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . = %12.2f\n' ...
    '     TIME FOR LOAD CASE SOLUTIONS  . . . . . . . . . . = %12.2f\n\n' ...
    '      T O T A L   S O L U T I O N   T I M E  . . . . . = %12.2f\n'], ...
    time(1), time(2), time(3), time(4));

fprintf(['\n' ...
    ' S O L U T I O N   T I M E   L O G   I N   S E C\n\n' ...
    '     TIME FOR INPUT PHASE  . . . . . . . . . . . . . . = %12.2f\n' ...
    '     TIME FOR CALCULATION OF STIFFNESS MATRIX  . . . . = %12.2f\n' ...
    '     TIME FOR LOAD CASE SOLUTIONS  . . . . . . . . . . = %12.2f\n\n' ...
    '      T O T A L   S O L U T I O N   T I M E  . . . . . = %12.2f\n'], ...
    time(1), time(2), time(3), time(4));

fclose(cdata.IIN);
fclose(cdata.IOUT);
end
end