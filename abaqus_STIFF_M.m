clc;clear;close all;
Dir = 'D:\SIMULIA\Temp\';
% Hed ='Job-C3D8-cantilever-beam-rho7800-dyna-matrix';
% Hed ='Job-C3D20-cantilever-beam-rho7800-matrix';
% Hed ='Job-C3D8-dyna-impicit-matrix';
Hed ='Job-C3D8-canbeam-rho7800-dyna-800s-matrix';

MASS_ID = [Hed,'_MASS2.mtx'];
STIFF_ID = [Hed,'_STIF2.mtx'];

fid_mass = fopen([Dir,MASS_ID],'r');
Mass = GetMatrix(fid_mass);

fid_stiff = fopen([Dir,STIFF_ID],'r');
Stiff = GetMatrix(fid_stiff);

% bc_num = [1,3,5,7];
bc_num = [9,6,3,8,5,2,7,4,1];
% num = [9,130,6,116,3,131,117,115,8,123,5,105,2,124,106,104,7,125,4,107,1];
num_freedom_id = [(bc_num-1).*3+1, (bc_num-1).*3+2, (bc_num-1).*3+3];
Stiff_total = Stiff;
Stiff_total(num_freedom_id,:) = [];
Stiff_total(:,num_freedom_id) = [];

Mass_total = Mass;
Mass_total(num_freedom_id,:) = [];
Mass_total(:,num_freedom_id) = [];

rank_Stiff_total = rank(Stiff_total);
rank_Mass_total = rank(Mass_total);

% q_load = zeros(24,1);%行号节点编号，值 载荷大小
% q_load((6-1)*3+2)=-20;

% Q_total = q_load;
% Q_total(num_freedom_id,:) = [];

function matrix = GetMatrix(fid)
a=textscan(fid,'%d %d %d %d %f','Delimiter',',');%
n = max(a{1})*3;
matrix = zeros(n,n);
for i = 1:length(a{1})
        matrix( 3*(a{1}(i)-1)+a{2}(i) , 3*(a{3}(i)-1)+a{4}(i) )=a{5}(i);
        matrix( 3*(a{3}(i)-1)+a{4}(i) , 3*(a{1}(i)-1)+a{2}(i) )=a{5}(i);
end

end
