%* *****************************************************************
%* - Basic data class of STAPMAT                                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Storing variables which control the running of STAPMAT      *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.20                *
%*                                                                 *
%* *****************************************************************
classdef ControlData
    properties
        NUMNP;         % Total number of nodal points
                       % = 0 : Program stop

        NPAR;          % Element group control data
                       %   NPAR(1) - Element type
                       %             1 : Truss element
                       %   NPAR(2) - Number of elements
                       %   NPAR(3) - Number of different sets of material
                       %             and cross-sectional constants
        NUMEG;         % Total number of element groups, > 0
        NUMEM;        % Total number of elements, >0
        NLCASE;        % Number of load case (>0)
        LL;            % Load case number
        NLOAD;         % The number of concentrated loads applied in this load case

        DLTYPE;        % Dynamic Load Type，动力学载荷类型；1-正弦（包含常数）
        DLPAR;         % Dynamic Load Parameters，动力学载荷集参数

        MODEX;         % Solution mode: 0 - data check only; 1 - execution; 2 - MATLAB稀疏矩阵左除求解
        DSTIME;        % 动力学求解时间：从0到DSTIME
        DSMETH;        % 动力学求解方法：1-精细积分；2-模态叠加；3-Verlet；4-广义alpha

        TIM;           % Timing information
        HED;           % Master heading information for usr in labeling the output

        IIN;           % file pointer used for input
        IOUT;        % file pointer used for output
        IDAT_ANIM;        % file pointer used for output for Tecplot
        IDAT_CURV;        % file pointer used for output for Tecplot
    end
end
