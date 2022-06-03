% % Called by GetStiff.m
% % 
% % - Call procedures:
% %        InitC3D8() % 定义每个单元上的节点数和每个节点的自由度
% %        ReadC3D8()
% % 

function C3D8Stiff(NUMEG_ID)

% Init variables of the element
InitC3D8();

% Read Material and Elements
ReadC3D8(NUMEG_ID);

fprintf('Solution phase ...\n\n');

% calculate addresses of diagonal elements
Addres();

% Data check Or Solve
global cdata;
if (cdata.MODEX == 0) 
    cdata.TIM(3,:) = clock;
    cdata.TIM(4,:) = clock;
    cdata.TIM(5,:) = clock;
    return; 
end

% Assemble structure stiffness matrix
Assemble();

end % end of Solid3DStiff


% ----------------------- Functions -----------------------------------

% Init parameters of C3D8 element
function InitC3D8()

global sdata;
sdata.NNODE = 8; % 一个单元上的节点数
sdata.NDOF = 3; % 每个节点的自由度

end % end of function InitC3D8()


% Assemble structure stiffness matrix
function Assemble()

global sdata;
global cdata;
sdata.STIFF = zeros(sdata.NWK, 1, 'double');
sdata.MASS = zeros(sdata.NWK, 1, 'double');
MASS_DIRCTION = sdata.MASS_DIRCTION;

ELEII = sdata.ELEII;
ID = sdata.ID;

NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
E = sdata.E; nu = sdata.nu; NGaussian = sdata.NGaussian; LM = sdata.LM;rho = sdata.rho;

for N = 1:NUME
    MTYPE = MATP(N);
    rhoE = rho(MTYPE);

    [ng,ksi,eta,zeta,weight] = Get3DGaussianIntInfo( NGaussian );

    E0 = E(MTYPE);
    nu0 = nu(MTYPE);
    D0 = E0*(1-nu0) / ( (1+nu0) * (1-2*nu0) );
    D1 = nu0/(1-nu0);
    D2 = (1-2*nu0) / ( 2*(1-nu0) );
    D = D0 * [1, D1, D1, 0, 0, 0;
              D1, 1, D1, 0, 0, 0;
              D1, D1, 1, 0, 0, 0;
               0, 0, 0, D2, 0, 0;
               0, 0, 0, 0, D2, 0;
               0, 0, 0, 0, 0, D2]; % sigma = D * epsilon

    node_coor = XYZ(:,N); % 该单元上的节点的XYZ坐标
    node_id = ELEII(N,:);
    dof_id = zeros(length(node_id)*3,1);
    for i_node_id = 1:length(node_id)
        dof_id(i_node_id*3-2) = ID(1,node_id(i_node_id));
        dof_id(i_node_id*3-1) = ID(2,node_id(i_node_id));
        dof_id(i_node_id*3)   = ID(3,node_id(i_node_id));
    end
%     dof_id

%     node_coor = zeros(sdata.NNODE*3,1); % 一个单元上的节点的xyz坐标
%     node_coor(1:3:end) = XYZ(1:3:end,N);
%     node_coor(2:3:end) = XYZ(2:3:end,N);
%     node_coor(3:3:end) = XYZ(3:3:end,N);

    Ke = zeros(sdata.NNODE*3,sdata.NNODE*3); % 单元刚度阵
    Me_x = zeros(sdata.NNODE*3,sdata.NNODE*3); % 单元协调质量阵
    dV = 0;
    for i = 1:ng
        for j = 1:ng
            for k = 1:ng
                [Shape,Jacobi,B] = C3D8NJB(node_coor,ksi(i),eta(j),zeta(k));
                Ke = Ke + weight(i)*weight(j)*weight(k) * (B') * D * B * det(Jacobi);
                Me_x = Me_x + weight(i)*weight(j)*weight(k) * rhoE * (Shape') * Shape;
                dV = dV + det(Jacobi);
            end
        end
    end
% % % % 直接组装协调质量阵    
%     for i_Mex = 1:length(dof_id)
%         for j_Mex = 1:length(dof_id)
%             if (dof_id(i_Mex) ~= 0) && (dof_id(j_Mex) ~= 0)
%                 MASS_DIRCTION( dof_id(i_Mex), dof_id(j_Mex) ) = MASS_DIRCTION( dof_id(i_Mex), dof_id(j_Mex) ) + Me_x( i_Mex, j_Mex );
%             end
%         end
%     end
% % % % 
    dV = mean(dV);
    alpha = rhoE * dV / sum(diag(Me_x));
    Me = zeros(sdata.NNODE*3,sdata.NNODE*3); % 单元集中质量阵
    for i = 1:sdata.NNODE*3
%         Me(i,i) = sum(Me_x(i,:));
        Me(i,i) = alpha * Me_x(i,i);
    end

%     Me = Me_x; %协调质量阵，注释掉则使用集中质量阵

%   SRC/Mechanics/ADDBAN.m
    ADDBAN(Ke,Me,LM(:,N));
end
sdata.MASS_DIRCTION = MASS_DIRCTION;

% The third time stamp
cdata.TIM(3, :) = clock;

end % end of function Assemble()
