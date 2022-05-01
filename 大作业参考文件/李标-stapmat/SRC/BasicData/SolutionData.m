%* *****************************************************************
%* - Basic data class of STAPMAT                                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Storing variables used in solving process                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.20                *
%*                                                                 *
%* *****************************************************************

classdef SolutionData
    properties (Constant)
        % Gauss coord, 1D to 3D
        GC1 = double(0.0);
        GC2 = double([1/3,-1/3]);
        GC3 = double([sqrt(0.6), 0.0, -sqrt(0.6)]);
        % Gauss weight, 1D to 3D
        GW1 = double(2.0);
        GW2 = double([1.0, 1.0]);
        GW3 = double([5.0/9.0, 8.0/9.0, 5.0/9.0]);
    end
    properties
        % Basic data
        ID;       % int, ID(3, NUMNP), Boundary condition codes (0=free, 1=deleted)
        IDOrigin; % int, backups of ID after computing of NEQ
        X;        % double, X(NUMNP), X coordinates
        Y;        % double, Y(NUMNP), Y coordinates
        Z;        % double, Z(NUMNP), Z coordinates
        R;        % double, R(NEQ), Load vector
        NOD;      % int, NOD(NLOAD), Node number to which this load is applied (1~NUMNP)
        IDIRN;    % int, IDIRN(NLOAD), Degree of freedom number for this load component
                  %                     1 : X-direction;
                  %                     2 : Y-direction;
                  %                     3 : Z-direction;
        FLOAD;    % double, FLOAD(NLOAD), Magnitude of load

        
        
        % Element data
        NUME;     % int, number of elements
        NNODE;    % int, 单元节点数
        NINIP;    % int, number of integration points in an element
        NDOF;     % int, 位移自由度
        NSTIFF;   % int, the number of number in element stiffness matrix
        XYZ;      % double, XYZ(3*NNODE, NUME), element position
        
        InitCoord;  % double array, integration coordinates
        InitWeight; % double array, integration weights
        
        % Material data
        NUMMAT;     % int, the number of types of material 
        E;          % double array, Young's Modulus
        density;    % 质量密度
        PR;         % double array, possion ratio
        AREA;       % double array, cross-sectional constants
        MATP;       % int, MATP(NUME), types of elements
        NINT        % int, GAUSS NUMERICAL INTEGRATION ORDER
        THIC        % double array,THICKNESS OF ELEMENT
        
        % Solve data
        NEQ;      % int, Number of equations
        NWK;      % Number of matrix elements
        MK;       % Maximum half bandwidth
        MHT;      % int, MHT(NEQ), 每一自由度在总刚度阵中到对角线距离
        LM;       % int, LM(6, NUME), Connectivity matrix
        MAXA;     % int, MAXA(NEQ),对角线元素的地址
        STIFF;    % double ,STIFF(NWK), store the elements of stiffness matrix
        MASS;    %协调质量阵
        CHNOD;    %载荷节点坐标
        N;        %形函数
        BN;       %形函数导数
        masschoose; %生成质量矩阵的方式，1和3，1-是协调质量阵行和，3是节点积分。
        radius_way; %轴对称单元的半径形式，1是高斯点的半径，2是平均半径

        % Result data
        DIS;      % double, DIS(NEQ, NLCASE), Displacement of nodes
        V;        % double, 节点速度
        STRAIN;   % double, STRAIN(NEQ, NLCASE), Strain
        STRESS;   % double, STRESS(NEQ, NLCASE), Stress
        maxMises; % double, 最大米塞斯应力和单元编号
    end
end