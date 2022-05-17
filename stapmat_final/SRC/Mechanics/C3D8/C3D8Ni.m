% % C3D8单元形函数N_i
function [N_i] = C3D8Ni(ksi, eta, zeta)
%C3D8Ni 此处显示有关此函数的摘要
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

