% % 矩形单元形函数导数 B = LN
function B = PlaneLShape(Jacobi,ksi,eta)
    dN1dksi = -1/4*(1-eta);
    dN1deta = -1/4*(1-ksi);
    dN2dksi = 1/4*(1-eta);
    dN2deta = -1/4*(1+ksi);
    dN3dksi = 1/4*(1+eta);
    dN3deta = 1/4*(1+ksi);
    dN4dksi = -1/4*(1+eta);
    dN4deta = 1/4*(1-ksi);
    B0 = [dN1dksi, dN2dksi, dN3dksi, dN4dksi;
          dN1deta, dN2deta, dN3deta, dN4deta];
    B1 = Jacobi\B0;
    B = zeros(3,8);
    B(1,1:2:8) = B1(1,:);
    B(2,2:2:8) = B1(2,:);
    B(3,1:2:8) = B1(2,:);
    B(3,2:2:8) = B1(1,:);
end