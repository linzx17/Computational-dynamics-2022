% % C3D20单元形函数N，Jacobi矩阵，应变矩阵B=L*N
function [N,Jacobi,B] = C3D20NJB(node_coor, ksi, eta, zeta)

ksi_i  = [-1, 1, 1, -1, -1, 1, 1, -1, 0, 1, 0, -1, 0, 1, 0, -1, -1, 1, 1, -1];
eta_i  = [-1, -1, 1, 1, -1, -1, 1, 1, -1, 0, 1, 0, -1, 0, 1, 0, -1, -1, 1, 1];
zeta_i = [-1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 0, 0, 0, 0];

ksi_0 = ksi_i * ksi;
eta_0 = eta_i * eta;
zeta_0 = zeta_i * zeta;

%% 形函数N
N_i = zeros(1,20);

N_i(1:8) = 1/8 * (1+ksi_0(1:8)) .* (1+eta_0(1:8)) .* (1+zeta_0(1:8)) .* (ksi_0(1:8) + eta_0(1:8) + zeta_0(1:8) - 2);

for i = [9, 11, 15, 13]
    N_i(i) = 1/4 * (1-ksi^2) * (1+eta_0(i)) * (1+zeta_0(i));
end

for i = [10, 12, 16, 14]
    N_i(i) = 1/4 * (1+ksi_0(i)) * (1-eta^2) * (1+zeta_0(i));
end

for i = [17, 18, 19, 20]
    N_i(i) = 1/4 * (1+ksi_0(i)) * (1+eta_0(i)) * (1-zeta^2);
end

N = zeros(3,60); % 形函数
N(1, 1:3:end) = N_i;
N(2, 2:3:end) = N_i;
N(3, 3:3:end) = N_i;

%% Jacobi 矩阵
dN = zeros(3,20);

dN(1,1:8) = 1/8 * ksi_i(1:8) .* (1+eta_0(1:8)) .* (1+zeta_0(1:8)) .* (ksi_0(1:8) + eta_0(1:8) + zeta_0(1:8) - 2) + ...
            1/8 * (1+ksi_0(1:8)) .* (1+eta_0(1:8)) .* (1+zeta_0(1:8)) .* ksi_i(1:8);
dN(2,1:8) = 1/8 * eta_i(1:8) .* (1+ksi_0(1:8)) .* (1+zeta_0(1:8)) .* (ksi_0(1:8) + eta_0(1:8) + zeta_0(1:8) - 2) + ...
            1/8 * (1+ksi_0(1:8)) .* (1+eta_0(1:8)) .* (1+zeta_0(1:8)) .* eta_i(1:8);
dN(3,1:8) = 1/8 * zeta_i(1:8) .* (1+ksi_0(1:8)) .* (1+eta_0(1:8)) .* (ksi_0(1:8) + eta_0(1:8) + zeta_0(1:8) - 2) + ...
            1/8 * (1+ksi_0(1:8)) .* (1+eta_0(1:8)) .* (1+zeta_0(1:8)) .* zeta_i(1:8);

for i = [9, 11, 15, 13]
    dN(1,i) = 1/4 * (-2*ksi) * (1+eta_0(i)) * (1+zeta_0(i));
    dN(2,i) = 1/4 * (1-ksi^2) * eta_i(i) * (1+zeta_0(i));
    dN(3,i) = 1/4 * (1-ksi^2) * (1+eta_0(i)) * zeta_i(i);
end

for i = [10, 12, 16, 14]
%     N_i(i) = 1/4 * (1+ksi_0(i)) * (1-eta^2) * (1+zeta_0(i));
    dN(1,i) = 1/4 * ksi_i(i) * (1-eta^2) * (1+zeta_0(i));
    dN(2,i) = 1/4 * (1+ksi_0(i)) * (-2*eta) * (1+zeta_0(i));
    dN(3,i) = 1/4 * (1+ksi_0(i)) * (1-eta^2) * zeta_i(i);
end

for i = [17, 18, 19, 20]
%     N_i(i) = 1/4 * (1+ksi_0(i)) * (1+eta_0(i)) * (1-zeta^2);
    dN(1,i) = 1/4 * ksi_i(i) * (1+eta_0(i)) * (1-zeta^2);
    dN(2,i) = 1/4 * (1+ksi_0(i)) * eta_i(i) * (1-zeta^2);
    dN(3,i) = 1/4 * (1+ksi_0(i)) * (1+eta_0(i)) * (-2*zeta);
end

coor = zeros(20,3); % 第1,2,3列分别表示节点的xyz坐标
coor(:,1) = node_coor(1:3:end);
coor(:,2) = node_coor(2:3:end);
coor(:,3) = node_coor(3:3:end);

Jacobi = dN * coor;

%% B = L * N
dNdxyz = Jacobi \ dN;
B = zeros(6,60);

[B(1,1:3:end), B(4,2:3:end), B(6,3:3:end)] = deal( dNdxyz(1,:) );
[B(2,2:3:end), B(4,1:3:end), B(5,3:3:end)] = deal( dNdxyz(2,:) );
[B(3,3:3:end), B(5,2:3:end), B(6,1:3:end)] = deal( dNdxyz(3,:) );

end