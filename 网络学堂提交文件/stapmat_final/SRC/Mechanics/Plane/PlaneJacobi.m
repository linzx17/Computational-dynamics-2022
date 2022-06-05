% % 根据四边形节点坐标计算jacobi矩阵
function Jacobi = PlaneJacobi(node_coor,ksi,eta)
        
    dN = 1/4 * [-(1-eta), 1-eta, 1+eta, -(1+eta);
                -(1-ksi), -(1+ksi), 1+ksi, 1-ksi];
    
    node_xy = [node_coor(1:2:7), node_coor(2:2:8)];

    Jacobi = dN * node_xy;
    
end