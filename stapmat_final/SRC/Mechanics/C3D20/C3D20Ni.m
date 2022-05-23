%-Called by :C3D20Stress.m

% % C3D20单元形函数N_i
function [N_i] = C3D20Ni(ksi, eta, zeta)
%下面(ksi_i(1),eta_i(1),zeta_i(1))...(ksi_i(20),eta_i(20),zeta_i(20)))表示的即为strain_n_1(1,:)中依次求出的节点顺序
ksi_i  = [-1, 1, 1, -1, -1, 1, 1, -1, 0, 1, 0, -1, 0, 1, 0, -1, -1, 1, 1, -1];
eta_i  = [-1, -1, 1, 1, -1, -1, 1, 1, -1, 0, 1, 0, -1, 0, 1, 0, -1, -1, 1, 1];
zeta_i = [-1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 0, 0, 0, 0];

ksi_0 = ksi_i * ksi;
eta_0 = eta_i * eta;
zeta_0 = zeta_i * zeta;

%形函数矩阵分量
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

end