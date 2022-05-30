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
        GC2 = double([-1/sqrt(3),1/sqrt(3)]);
        GC3 = double([sqrt(0.6), 0.0, -sqrt(0.6)]);
        % Gauss weight, 1D to 3D
        GW1 = double(2.0);
        GW2 = double([1.0, 1.0]);
        GW3 = double([5.0/9.0, 8.0/9.0, 5.0/9.0]);
        % Rayleigh阻尼系数：阻尼矩阵 C = Calpha*M + Cbeta*K;
        Calpha = 0; 
        Cbeta = 0;
    end
    properties
        NGaussian;
        GC;
        GW;
        % Basic data
        ID;       % int, ID(3, NUMNP), Boundary condition codes (0=free, 1=deleted), ID数组的第1、2、3行分别表示xyz方向自由度，第j(1<=j<=cdata.NUMNP)列表示第j个节点
        IDOrigin; % int, backups of ID after computing of NEQ
        X;        % double, X(NUMNP), X coordinates
        Y;        % double, Y(NUMNP), Y coordinates
        Z;        % double, Z(NUMNP), Z coordinates
        ELEII;   %int, ELEII(NUMEM,NNODE), element information
        NUMEGEM;   % int, NUMEGEM(NUMEG),Numid of elements in one element group;
        NODE_FLAG; % double, NODE_FLAG(NUMNP,NLCASE),In all element, the times of the node appears

        R;        % double, R(NEQ), Load vector
        NOD;      % int, NOD(NLOAD), Node number to which this load is applied (1~NUMNP)
        IDIRN;    % int, IDIRN(NLOAD), Degree of freedom number for this load component
                  %                     1 : X-direction;
                  %                     2 : Y-direction;
                  %                     3 : Z-direction;
        FLOAD;    % double, FLOAD(NLOAD), Magnitude of load

        
        
        % Element data
        NUME;     % int, number of elements in this element group
        NNODE;    % int, number of nodes in an element
        NINIP;    % int, number of integration points in an element
        NDOF;     % int, the DOF of displacement
        NSTIFF;   % int, the number of number in element stiffness matrix
        XYZ;      % double, XYZ(3*NNODE, NUME), element position
        
        InitCoord;  % double array, integration coordinates
        InitWeight; % double array, integration weights
        
        % Material data
        NUMMAT;     % int, the number of types of material 
        E;          % double array, Young's Modulus
        nu;         % double array, possion ratio
        rho;        % double array, density
        AREA;       % double array, cross-sectional constants
        MATP;       % int, MATP(NUME), types of elements
        
        % Solve data
        NEQ;      % int, Number of equations
        NWK;      % Number of matrix elements
        MK;       % Maximum half bandwidth
        MHT;      % int, MHT(NEQ), Vector of column heights
        LM;       % int, LM(6, NUME), Connectivity matrix
        MAXA;     % int, MAXA(NEQ)
        STIFF;    % double ,STIFF(NWK), store the elements of stiffness matrix
        MASS;     % double ,MASS(NWK), store the elements of mass matrix
        SPSTIFF;  % 刚度阵的matlab稀疏矩阵
        SPMASS;   % 质量阵的matlab稀疏矩阵
        FRE_EIG;  % 调用eig函数求得的所有固有频率(Hz)
        FREQUENCY;% 固有频率(Hz) = sqrt(特征值)/2/pi
        EIGVALUE; % 特征值
        EIGVECTOR;% 阵型
        EIGDISP;  % 阵型对应的节点位移
        
        % Result data
        DIS;      % double, DIS(NEQ, NLCASE), Displacement of nodes
        DISP;    % double, DIS(NUMNP, NLCASE), Displacement of nodes
        STRAIN;   % double, STRAIN(NUMNP,6, NLCASE), Strain
        STRESS;   % double, STRESS(NUMNP,6, NLCASE), Stress
        DYNADIS;    % 时间积分计算结果：每一列为一个自由度，行为时间
        DYNADT;     % DYDIS结果中对应的时间序列
        
    end
end