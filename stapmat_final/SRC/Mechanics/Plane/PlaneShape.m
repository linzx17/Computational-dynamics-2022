% % 矩形单元形函数
function N = PlaneShape(ksi,eta)
    N1 = 1/4*(1-ksi)*(1-eta);
    N2 = 1/4*(1+ksi)*(1-eta);
    N3 = 1/4*(1+ksi)*(1+eta);
    N4 = 1/4*(1-ksi)*(1+eta);
    N = zeros(2,8);
    N(1:2:8,1) = [N1,N2,N3,N4];
    N(2:2:8,2) = [N1,N2,N3,N4];
end