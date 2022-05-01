% % Call procedures:
% % Called by: SRC/Solver/GetStress.m
% % 在out文件中输出高斯积分点上的应力应变

function PlaneStress(NUM,NG)


% Get global data
global cdata;
global sdata;

IOUT = cdata.IOUT;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
E = sdata.E; nu = sdata.nu; LM = sdata.LM;
U = sdata.DIS(:, NUM);

fprintf(IOUT, ['\n\n  S T R A I N AND S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '  ELEMENT        GAUSS POINT                                    STRAIN                                                                       STRESS\n' ...
    '  NUMBER            NUMBER            EPSILON_X         EPSILON_Y        GAMMA_XY              SIGMA_X            SIGMA_Y             TAU_XY\n'], NG);

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

    ae = zeros(8,1);% 节点位移
    for i = 1:8
        if LM(i,N) > 0
            ae(i) = U(LM(i,N));
        end
    end
    ksi = sdata.GC;
    eta = sdata.GC;
    weight = sdata.GW;
    ng = sdata.NG;
    strain_e = zeros(3,ng*ng);%高斯积分点上的应变
    stress_e = zeros(3,ng*ng);%高斯积分点上的应力

    for i = 1:ng
        for j = 1:ng
            Jacobi = PlaneJacobi(node_coor,ksi(j),eta(i));
            B = PlaneLShape(Jacobi,ksi(j),eta(i));
            gpnum = (i-1)*ng+j;%高斯点计数
            strain_e(:,gpnum) = B*ae;
            stress_e(:,gpnum) = D*strain_e(:,gpnum);
            fprintf(IOUT,'       %d                        %d                %13.6e    %13.6e    %13.6e        %13.6e    %13.6e    %13.6e\n', ...
                N, gpnum, strain_e(1,gpnum), strain_e(2,gpnum), strain_e(3,gpnum), stress_e(1,gpnum), stress_e(2,gpnum), stress_e(3,gpnum));
        end
    end

end


end