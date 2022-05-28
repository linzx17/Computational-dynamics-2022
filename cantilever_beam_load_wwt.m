clc;clear;close all;
%SI(mm)
%mm N tonne s MPa mJ tonne/mm^3
%悬臂梁纯弯曲应力状态，需要在梁端截面施加沿高度方向线性分布的面载荷[N/mm^2]
%弯曲力矩为-100000[N·mm]
M = 100000;%[N·mm]
S = 5*10;%[mm^2]
L = 5;%[mm]
q_surfload = M/S/L;%[N/mm^2]
%对于C3D20单元而言，均匀面载荷为拉力，大小为1N
%等价于对各个角节点施加-1/12N的力，对边中节点施加1/3N的力
%对于一个单元而言，边长为a
a = 5;%[mm]
q_node = a^2*q_surfload*(-1/12);
q_egde = a^2*q_surfload*(1/3);

q_load = zeros(21,1);%行号节点编号，值 载荷大小

num = [1 3 9 11];
q_load(num,1) = q_load(num,1) + q_node.*ones(4,1);
num = [3 5 13 11];
q_load(num,1) = q_load(num,1) + q_node.*ones(4,1);
num = [9 11 19 17];
q_load(num,1) = q_load(num,1) - q_node.*ones(4,1);
num = [11 13 21 19];
q_load(num,1) = q_load(num,1) - q_node.*ones(4,1);

num = [2 7 10 6];
q_load(num,1) = q_load(num,1) + q_egde.*ones(4,1);
num = [4 8 12 7];
q_load(num,1) = q_load(num,1) + q_egde.*ones(4,1);
num = [10 15 18 14];
q_load(num,1) = q_load(num,1) - q_egde.*ones(4,1);
num = [12 16 20 15];
q_load(num,1) = q_load(num,1) - q_egde.*ones(4,1);

bc_node= [9,130,6,116,3,131,117,115,8,123,5,105,2,124,106,104,7,125,4,107,1];
load_node = [93 310 96 320 99 311 309 319 92 303 95 316 98 304 302 315 91 301 94 314 97];
Q = zeros(963,1);
for i = 1: 21
    Q((load_node(i)-1)*3+1,1) = q_load(i,1);
end
num_freedom_id = [(bc_node-1).*3+1, (bc_node-1).*3+2, (bc_node-1).*3+3];
Q_total = Q;
Q_total(num_freedom_id,:) = [];