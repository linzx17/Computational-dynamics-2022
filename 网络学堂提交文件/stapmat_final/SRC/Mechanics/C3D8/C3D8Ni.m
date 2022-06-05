%-Called by :C3D8Stress.m

% % C3D8单元形函数N_i
function [N_i] = C3D8Ni(ksi, eta, zeta)
%下面(ksi_i(1),eta_i(1),zeta_i(1))...(ksi_i(8),eta_i(8),zeta_i(8)))表示的即为strain_n_1(1,:)中依次求出的节点顺序
ksi_i = [-1; 1; 1; -1; -1; 1; 1; -1];
eta_i = [-1; -1; 1; 1; -1; -1; 1; 1];
zeta_i = [-1; -1; -1; -1; 1; 1; 1; 1];

ksi_0 = ksi_i * ksi;
eta_0 = eta_i * eta;
zeta_0 = zeta_i * zeta;

%形函数矩阵分量
N_i = zeros(1,8);
for i = 1:8
    N_i(i) = 1/8 * (1+ksi_0(i)) * (1+eta_0(i)) * (1+zeta_0(i));
end

end

