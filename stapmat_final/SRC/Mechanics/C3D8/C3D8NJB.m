% % C3D8单元形函数N，Jacobi矩阵，应变矩阵B=L*N
function [N,Jacobi,B] = C3D8NJB(node_coor, ksi, eta, zeta)

ksi_i = [-1; 1; 1; -1; -1; 1; 1; -1];
eta_i = [-1; -1; 1; 1; -1; -1; 1; 1];
zeta_i = [-1; -1; -1; -1; 1; 1; 1; 1];

% ksi_i = [1; 1; 1; 1; -1; -1; -1; -1];
% eta_i = [-1; 1; 1; -1; -1; 1; 1; -1];
% zeta_i = [1; 1; -1; -1; 1; 1; -1; -1];

ksi_0 = ksi_i * ksi;
eta_0 = eta_i * eta;
zeta_0 = zeta_i * zeta;

%% 形函数
N_i = zeros(1,8);
for i = 1:8
    N_i(i) = 1/8 * (1+ksi_0(i)) * (1+eta_0(i)) * (1+zeta_0(i));
end

N = zeros(3,8*3); % 形函数
N(1,1:3:end) = N_i;
N(2,2:3:end) = N_i;
N(3,3:3:end) = N_i;

%% Jacobi矩阵
dN = zeros(3,8);
for i = 1:8
    dN(1,i) = ksi_i(i) * (1+eta_0(i)) * (1+zeta_0(i)) / 8; % dNdksi
    dN(2,i) = eta_i(i) * (1+ksi_0(i)) * (1+zeta_0(i)) / 8; % dNdeta
    dN(3,i) = zeta_i(i) * (1+ksi_0(i)) * (1+eta_0(i)) / 8; % dNdzeta
end

coor = zeros(8,3); % 第1,2,3列分别表示节点的xyz坐标
coor(:,1) = node_coor(1:3:end);
coor(:,2) = node_coor(2:3:end);
coor(:,3) = node_coor(3:3:end);

Jacobi = dN * coor; % 

%% B = L*N
dNdxyz = Jacobi \ dN;
B = zeros(6,24);

[B(1,1:3:end), B(4,2:3:end), B(6,3:3:end)] = deal( dNdxyz(1,:) );
[B(2,2:3:end), B(4,1:3:end), B(5,3:3:end)] = deal( dNdxyz(2,:) );
[B(3,3:3:end), B(5,2:3:end), B(6,1:3:end)] = deal( dNdxyz(3,:) );

end