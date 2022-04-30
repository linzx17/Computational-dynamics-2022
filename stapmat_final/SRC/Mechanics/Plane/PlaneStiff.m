
function PlaneStiff()

% Init variables of the element
InitPlane();

% Read Material and Elements
ReadPlane();

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




end

% ----------------------- Functions -----------------------------------

% Init parameters of truss element
function InitPlane()
global sdata;
sdata.NNODE = 4;% 一个单元上的节点数
sdata.NDOF = 2;% 每个节点的自由度

end

% Assemble structure stiffness matrix
function Assemble()
global sdata;
global cdata;
sdata.STIFF = zeros(sdata.NWK, 1, 'double');

NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
E = sdata.E; nu = sdata.nu; LM = sdata.LM;
if sdata.NG == 1
    sdata.GC = sdata.GC1;
    sdata.GW = sdata.GW1;
elseif sdata.NG == 2
    sdata.GC = sdata.GC2;
    sdata.GW = sdata.GW2;
else
    sdata.GC = sdata.GC3;
    sdata.GW = sdata.GW3;
end
ng = sdata.NG;
ksi = sdata.GC;
eta = sdata.GC;
weight = sdata.GW;
for N = 1:NUME
    MTYPE = MATP(N);
    E0 = E(MTYPE);
    nu0 = nu(MTYPE);
    D0 = E0/(1-nu0^2);
    D = D0*[1, nu0, 0;
        nu0, 1, 0;
        0, 0, (1-nu0)/2];% sigma = D * epsilon

    node_coor = zeros(8,1);% 四个节点的xy坐标
    node_coor(1:2:end) = XYZ(1:3:end,N);
    node_coor(2:2:end) = XYZ(2:3:end,N);

    Ke = zeros(8,8);% 单元刚度阵
    for i = 1:ng
        for j = 1:ng
            Jacobi = PlaneJacobi(node_coor,ksi(i),eta(j));
            B = PlaneLShape(Jacobi,ksi(i),eta(j));
            Ke = Ke + weight(i)*weight(j)*B'*D*B*det(Jacobi);
        end
    end

%   SRC/Mechanics/ADDBAN.m
    ADDBAN(Ke,LM(:,N));
end

% The third time stamp
cdata.TIM(3, :) = clock;

end
