%%将abaqus输入文件更换为matlab输入文件
%提取abaqus的节点和单元信息
%abaqus节点信息 包括 编号 xyz坐标 四维数组
clc;clear;close all;
Dir = 'D:\SIMULIA\Temp';

hed = 'C3D20_cantilever_beam_matlab';                  %文件头
hed_Flag = 2;                                                            % 1-C3D8  2-C3D20
numnp = 321;                                                             %节点个数
numel = 40;                                                                 %单元个数
ID(1:3,numnp)=0;
ID(1:3,[9,130,6,116,3,131,117,115,8,123,5,105,2,124,106,104,7,125,4,107,1])=1; 
%添加约束

abaqus_inp_ID = ['\Job_',hed,'.inp'];%abaqus路径
matlab_in_ID = ['\Job_',hed,'.in'];%matlab路径

fileID_abaqus= fopen([Dir,abaqus_inp_ID]);
node_data=textscan(fileID_abaqus,'%d %f %f %f',numnp,'Delimiter',',','headerlines',9);
fclose(fileID_abaqus);
fileID_abaqus= fopen([Dir,abaqus_inp_ID]);
if hed_Flag == 1
    element_data=textscan(fileID_abaqus,'%d %d %d %d %d %d %d %d %d',numel,'Delimiter',',','headerlines',10+numnp);
elseif hed_Flag == 2
    element_data=textscan(fileID_abaqus,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',numel,'Delimiter',',','headerlines',10+numnp);
end
fclose(fileID_abaqus);

numnp = length(node_data{1});%节点个数
numel = length(element_data{1});%单元个数

numeg = 1;%单元组总数
numem = 1;%单元总数
%length(npar(:,1)) = numeg
if hed_Flag == 1
    npar(1,1) = 3;                                                              %单元类型 3-C3D8 4-C3D20
elseif hed_Flag == 2
    npar(1,1) = 4;
else
    disp('单元类型npar1检查');
end
npar(1,2) = numel;%本单元组中单元总数
npar(1,3) = 1;%本单元组中不同材料/截面性质组数
    %length(material(:,1)) = npar(1,3)
    material(1,1) = 210000;                                              %E(1)   杨氏模量
    material(1,2) = 0.3;                                                     %MU(1)  泊松比
    material(1,3) = 7.8e-9;                                                   %RHO(1) 密度 1*10^3[kg/(mm^3)]
npar(1,4) = 2;                                                                  %高斯积分阶次

nlcase = 1;%载荷工况数
    nload(1) = 16;%集中载荷的个数
    %length(load_data(:,1)) = nload(1)
    load_data(1:nload(1),2) = [93 99 310 320 311 319 96 309 304 315 301 314 302 91 97 94];
    %集中载荷作用的节点号
    load_data(1:nload(1),3) = 1.*ones(16,1);                                       %载荷作用方向
    load_data(1:nload(1),4) = [-833.33 -833.33 3333.3 3333.3 3333.3 3333.3 ...
        -1666.7 6666.7 -3333.3 -3333.3 -3333.3 -3333.3 -6666.7 833.33 833.33 1666.7];                                   %载荷作用值
modex = 2;%求解模式


fileID_matlab = fopen([Dir,matlab_in_ID], 'w+');
%输出标题行
fprintf(fileID_matlab,hed);
fprintf(fileID_matlab,'\n');
%输出控制行
fprintf(fileID_matlab,' %4d %4d %4d %4d %4d\n', numnp,numeg,numem,nlcase,modex);
%输出节点数据
for i = 1:numnp
    fprintf(fileID_matlab,' %4d %4d %4d %4d  %17.9E %17.9E %17.9E\n',i,ID(1,i),ID(2,i),ID(3,i),node_data{2}(i),node_data{3}(i),node_data{4}(i));
end
%输出载荷数据
for i = 1:nlcase
    fprintf(fileID_matlab,' %4d %4d\n',i,nload(i));
    for j = 1:nload(i)
        fprintf(fileID_matlab,' %4d %4d %9d\n',load_data(j,2),load_data(j,3),load_data(j,4));
    end
end
%三维单元数据
for i = 1:numeg
    %单元组控制数据
    fprintf(fileID_matlab,' %4d %4d %4d %4d\n',npar(i,1),npar(i,2),npar(i,3),npar(i,4));
    %材料/截面性质数据
    for j = 1:npar(i,3)
        fprintf(fileID_matlab,' %4d %9d %4.1f %9.1E\n',j,material(j,1),material(j,2),material(j,3));
    end
    %单元数据
    for k = 1:numel
        if length(element_data)-1 == 8 %C3D8
            fprintf(fileID_matlab,' %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n',k...
                ,element_data{2}(i),element_data{3}(i),element_data{4}(i),element_data{5}(i)...
                ,element_data{6}(i),element_data{7}(i),element_data{8}(i),element_data{9}(i)...
                ,numeg...
                );
        elseif length(element_data)-1 == 20 %C3D20
            fprintf(fileID_matlab,' %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n',k...
                ,element_data{2}(i),element_data{3}(i),element_data{4}(i),element_data{5}(i)...
                ,element_data{6}(i),element_data{7}(i),element_data{8}(i),element_data{9}(i)...
                ,element_data{10}(i),element_data{11}(i),element_data{12}(i),element_data{13}(i)...
                ,element_data{14}(i),element_data{15}(i),element_data{16}(i),element_data{17}(i)...
                ,element_data{18}(i),element_data{19}(i),element_data{20}(i),element_data{21}(i)...
                ,numeg...
                );
        else
            fprintf(fileID_matlab,'Unknown solid element type!\n');
            fprintf(fileID_matlab,'Check element_data!\n');
        end
    end
end
fprintf(fileID_matlab,'stop'); 