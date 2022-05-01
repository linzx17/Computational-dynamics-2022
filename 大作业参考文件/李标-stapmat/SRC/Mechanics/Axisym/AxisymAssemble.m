%* *****************************************************************
%* - Function of STAPMAT in stiffness phase                        *
%*                                                                 *
%* - Purpose:                                                      *
%*     Assemble structure stiffness matrix                         *
%*                                                                 *
%* - Call procedures:                                              *
%*     SRC/Mechanics/ADDBAN.m - ADDBAN()                           *
%*                                                                 *
%* - Called by :                                                   *
%*     SRC/Mechanics/StrainStiff.m                                 *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.21                *
%*                                                                 *
%* *****************************************************************
function AxisymAssemble()
global sdata;


sdata.STIFF = zeros(sdata.NDOF^2 * sdata.NNODE^2 * sdata.NUME, 3);
sdata.MASS = zeros(sdata.NDOF^2 * sdata.NNODE^2 * sdata.NUME, 3);
NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
EN = sdata.E; LM = sdata.LM; NINT = sdata.NINT;
PRN = sdata.PR; density = sdata.density;
masschoose = sdata.masschoose;
radius_way = sdata.radius_way;
    XX = zeros(2, 8, 'double');    %XX是坐标，从2*4改成2*8
    XG=[0.0  -0.5773502691896  -0.7745966692415  -0.8611363115941;...
        0.0   0.5773502691896   0.0              -0.3399810435849;... %高斯积分点的位置，每一列对应着不同的积分方案，我们轴对称单元用2*2的减缩积分就行了
        0.0   0.0               0.7745966692415   0.3399810435849;...
        0.0   0.0               0.0               0.8611363115941];
    WGT=[2.0, 1.0, 0.5555555555555, 0.3478548451375;...%高斯积分点的权重
        0.0, 1.0, 0.8888888888888, 0.6521451548625;...
        0.0, 0.0, 0.5555555555555, 0.6521451548625;...
        0.0, 0.0, 0.0,             0.3478548451375];

  intiN(NINT); %初始化N,BN
for N = 1:NUME
    S = zeros(16, 16, 'double');  %S是刚度矩阵，将8*8改成16*16
    M = zeros(16, 16, 'double');  %M是刚度矩阵，将8*8改成16*16

    MTYPE = MATP(N);
    NU = PRN(MTYPE);
    E = EN(MTYPE);
    for i=1:2
        for j=1:8   %将4改成了8，因为是8个节点
            XX(i,j) = XYZ((j-1)*2+i,N);  %xx是物理坐标，第一行是坐标，第二行是y坐标
        end
    end   %检查了下XX没问题
    
    %生成单元的平均半径
    sumrow = sum(XX,2); %中介，为了计算平均半径利用一下
    rave = sumrow(1,1)/8;  %轴对称单元比较特殊，需要单元的平均半径,这里忘记除8了，但是不影响秩

    D=E/(1+NU)/(1-2*NU)*[1-NU NU 0 NU;...
                         NU 1-NU 0 NU;...
                         0 0 (1-2*NU)/2 0;...
                         NU NU 0 1-NU];                
    
    for LX=1:NINT  
     
        for LY=1:NINT
            
            [B,DET,rphy]=generateBAD(LX,LY,XX,radius_way);  %通过输入高斯积分点的位置，得到每个高斯积分点的B和行列式，这里的XX是有东西额，但是B没东西
            
            WT = WGT(LX,NINT)*WGT(LY,NINT)*DET*2*pi*rphy;%得到总权重
            S=S+B'*D*B*WT;    
        end   
    end
    
   %质量矩阵，要考虑积分点的权重
    for LX=1:3  
        RI = XG(LX,3); %将高斯积分点取出，这里的规律是先列扫描
     
        for LY=1:3
            SI=XG(LY,3);
            
            [ Nmatrix, DET] = generateN( RI, SI, XX ) ;  %通过输入高斯积分点的位置，得到每个高斯积分点的B和行列式，这里的XX是有东西额，但是B没东西
            WT = WGT(LX,3)*WGT(LY,3);
            M  = M + Nmatrix'* Nmatrix * DET * WT;%得到总权重
        end   
    end
    
   if masschoose == 1
    for LX=1:NINT  
        RI = XG(LX,NINT); %将高斯积分点取出，这里的规律是先列扫描
     
        for LY=1:NINT
            SI=XG(LY,NINT);
            
            [ Nmatrix, DET] = generateN( RI, SI, XX ) ;  %通过输入高斯积分点的位置，得到每个高斯积分点的B和行列式，这里的XX是有东西额，但是B没东西
            WT = WGT(LX,NINT)*WGT(LY,NINT);
            M  = M + Nmatrix'* Nmatrix * DET * WT;%得到总权重
        end   
    end
     M = 2 * pi * rave * density * M;
     M = diag(sum(M, 2)); 
        
        
    
    
    
   elseif masschoose == 2;
    M = 2 * pi * rave * density * M;
    midm1 = sum(M(:));
    M = diag(diag(M)); 
    M = M/(sum(M(:))/midm1);
    
    
   elseif masschoose == 3;
    for L = 1 : 8  
          [ Nmatrix, DET] = generateN( point_position(L, 1), point_position(L, 2), XX ) ;
          M  = M + Nmatrix'* Nmatrix * DET * weight_point(L);
    end
    M = 2 * pi * rave * density * M;
       
    end
%     
    %%这里是轴对称单元，所以单元的K矩阵需要乘以2*pi*r在循环中乘过了，以后可以乘以2pi，来减少程序的循环量
    ADDBANMASS(M, LM(:, N),N);
    ADDBAN(S, LM(:, N),N);

end    

% spy(sdata.STIFF)
end


function intiN(NINT)
global sdata;
sdata.N = zeros(8,NINT,NINT);
sdata.BN = zeros(2,8,NINT,NINT);
XG=[0.0  -0.5773502691896  -0.7745966692415  -0.8611363115941;...
    0.0   0.5773502691896   0.0              -0.3399810435849;...
    0.0   0.0               0.7745966692415   0.3399810435849;...
    0.0   0.0               0.0               0.8611363115941];
for LX=1:NINT
    RI = XG(LX,NINT); %将高斯积分点取出，这里的规律是先列扫描
    for LY=1:NINT
        SI=XG(LY,NINT);
        
        N(1) = 0.25*(1 - RI)*(1 - SI);
        N(2) = 0.25*(1 + RI)*(1 - SI);
        N(3) = 0.25*(1 + RI)*(1 + SI);
        N(4) = 0.25*(1 - RI)*(1 + SI);
        N(5) = 0.5*(1 - RI*RI)*(1 - SI); %将CQ4增添成了CQ8
        N(6) = 0.5*(1 + RI)*(1 - SI*SI);
        N(7) = 0.5*(1 - RI*RI)*(1 + SI);
        N(8) = 0.5*(1 - RI)*(1 - SI*SI);
        N(1) = N(1) - 0.5*N(5) - 0.5*N(8);%对原来的4个形函数进行修正
        N(2) = N(2) - 0.5*N(6) - 0.5*N(5);
        N(3) = N(3) - 0.5*N(7) - 0.5*N(6);
        N(4) = N(4) - 0.5*N(8) - 0.5*N(7);  %检查了N没有问题，一部分抽查
        
        N11 = 0.25*(1-SI)*(2*RI+SI);%1到4号节点是，左下角开始逆时针循环，5-8号是从底边的中点开始逆时针循环
        N12 = 0.25*(1-SI)*(2*RI-SI);
        N13 = 0.25*(1+SI)*(2*RI+SI);
        N14 = 0.25*(1+SI)*(2*RI-SI);
        
        N15 = -RI*(1-SI);  %增添了其余4个节点对psi求导
        N16 = 0.5*(1 - SI*SI);
        N17 = -RI*(1+SI);
        N18 = -0.5*(1 - SI*SI);
        
        N21 = 0.25*(1-RI)*(RI+2*SI);
        N22 = 0.25*(1+RI)*(-RI+2*SI);
        N23 = 0.25*(1+RI)*(RI+2*SI);
        N24 = 0.25*(1-RI)*(-RI+2*SI);
        
        N25 = -0.5*(1-RI*RI);  %增添了其余4个节点对eta求导
        N26 = -SI*(1+RI);
        N27 = 0.5*(1-RI*RI);
        N28 = -SI*(1-RI);
        %生成等参元的形函数的导数
        Bpsi = [N11, N12, N13, N14, N15, N16, N17, N18;...
            N21, N22, N23, N24, N25, N26, N27, N28];
        sdata.N(:,LX,LY) = N;
        sdata.BN(:,:,LX,LY) = Bpsi;
    end
end
end
