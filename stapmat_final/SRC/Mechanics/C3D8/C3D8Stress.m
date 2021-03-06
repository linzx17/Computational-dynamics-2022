% % Call procedures: C3D8NJB.m
% % Called by: SRC/Solver/GetStress.m
% % 在out文件中输出高斯积分点上的应力应变

% % 输入 载荷工况数NUM，单元组数NG
function C3D8Stress(NUM,NG)

% Get global data
global cdata;
global sdata;

IOUT = cdata.IOUT;
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ;
E = sdata.E; nu = sdata.nu; NGaussian = sdata.NGaussian; LM = sdata.LM;
U = sdata.DIS(:, NUM);

fprintf(IOUT, ['\n\n  S T R A I N AND S T R E S S  C A L C U L A T I O N S  F O R  ' ...
    'E L E M E N T  G R O U P %4d\n\n' ... 
    '  ELEMENT        GAUSS POINT                      \t\t\t     STRAIN    			            ' ...
    '                                                                                                        STRESS\n' ...
    '  NUMBER            NUMBER                  E11               E22               E33               E12' ...
    '               E23               E13               S11               S22               S33               S12               S23               S13\n'], NG);

for N = 1:NUME
    MTYPE = MATP(N);

    [ng,ksi,eta,zeta,~] = Get3DGaussianIntInfo( NGaussian );

    E0 = E(MTYPE);
    nu0 = nu(MTYPE);
    D0 = E0*(1-nu0) / ( (1+nu0) * (1-2*nu0) );
    D1 = nu0/(1-nu0);
    D2 = (1-2*nu0) / ( 2*(1-nu0) );
    D = D0 * [1, D1, D1, 0, 0, 0;
              D1, 1, D1, 0, 0, 0;
              D1, D1, 1, 0, 0, 0;
               0, 0, 0, D2, 0, 0;
               0, 0, 0, 0, D2, 0;
               0, 0, 0, 0, 0, D2]; % sigma = D * epsilon

    node_coor = zeros(sdata.NNODE*3,1); % 一个单元上的节点的xyz坐标
    node_coor(1:3:end) = XYZ(1:3:end,N);
    node_coor(2:3:end) = XYZ(2:3:end,N);
    node_coor(3:3:end) = XYZ(3:3:end,N);
    
    n_node_dof = sdata.NNODE*sdata.NDOF; % 一个单元上的节点自由度数
    ae = zeros(n_node_dof,1); % 第N个单元上的节点位移
    for i = 1:n_node_dof
        if LM(i,N) > 0
            ae(i) = U(LM(i,N));
        end
    end
    %单元中高斯积分点应力应变矩阵
    strain_e = zeros(6,ng^3);
    stress_e = zeros(6,ng^3);

    gpnum = 0; % 高斯积分点计数

    smooth_Ni = zeros(ng^3,8,'double');
    smooth_num = int8(0);

    for i = 1:ng
        for j = 1:ng
            for k = 1:ng
                [~,~,B] = C3D8NJB(node_coor,ksi(i),eta(j),zeta(k));
                smooth_num = smooth_num+1;
                smooth_Ni(smooth_num,:) = C3D8Ni(ksi(i),eta(j),zeta(k));
                gpnum = gpnum + 1;
                strain_e(:,gpnum) = B*ae.*[1 1 1 0.5 0.5 0.5]';
                stress_e(:,gpnum) = D*strain_e(:,gpnum);
                fprintf(IOUT,['    %d               %d             %13.6e     %13.6e     %13.6e     %13.6e     %13.6e     %13.6e' ...
                    '     %13.6e     %13.6e     %13.6e     %13.6e     %13.6e     %13.6e\n'], ...
                N, gpnum, strain_e(1,gpnum), strain_e(2,gpnum), strain_e(3,gpnum), strain_e(4,gpnum), strain_e(5,gpnum), strain_e(6,gpnum), ...
                stress_e(1,gpnum), stress_e(2,gpnum), stress_e(3,gpnum), stress_e(4,gpnum), stress_e(5,gpnum), stress_e(6,gpnum));
            end
        end
    end

    strain_n = (smooth_Ni \ strain_e')';
    stress_n = (smooth_Ni \ stress_e')';

    %strain_n(6,8) 其中6表示有6中应变，8表示有8个节点
    %下面ELEII的顺序依次表示节点顺序

    if (NG ~=1)
        element_id = sdata.NUMEGEM(NG-1)+N;
    else
        element_id = N;
    end

    node_id = sdata.ELEII(element_id,:);
%     sdata.NODE_FLAG(node_id,NUM) = sdata.NODE_FLAG(node_id,NUM)+ones(8,1,'int8');
    sdata.NODE_FLAG(node_id,NUM) = sdata.NODE_FLAG(node_id,NUM)+ones(8,1);

    sdata.STRAIN( sdata.ELEII(element_id,:) , : , NUM ) = sdata.STRAIN(sdata.ELEII(element_id,:)  , : , NUM)+strain_n(:,:)';
    sdata.STRESS( sdata.ELEII(element_id,:) , : , NUM ) = sdata.STRESS(sdata.ELEII(element_id,:)  , : , NUM)+stress_n(:,:)';

end