% % 矩形单元形函数
function N = PlaneShape(ksi,eta)
    N = zeros(2,8);
    N(1,1) = 1/4*(1-ksi)*(1-eta);
    N(1:3) = 1/4*(1+ksi)*(1-eta);
    N(1:5) = 1/4*(1+ksi)*(1+eta);
    N(1,7) = 1/4*(1-ksi)*(1+eta);
    N(2,2:2:8) = N(1,1:2:7);
end