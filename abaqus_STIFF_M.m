clc;clear;close all;
Dir = 'D:\SIMULIA\Temp\';
% Hed ='Job-C3D20-cantilever-beam-martix';
Hed ='Job-C3D20-cantilever-beam-rho7800-matrix';


MASS_ID = [Hed,'_MASS2.mtx'];
STIFF_ID = [Hed,'_STIF2.mtx'];

fid_mass = fopen([Dir,MASS_ID],'r');
Mass = GetMatrix(fid_mass);

fid_stiff = fopen([Dir,STIFF_ID],'r');
[Stiff,a_out] = GetMatrix(fid_stiff);

% rank_Mass = rank(Mass);
% rank_Stiff = rank(Stiff);

num = [9,130,6,116,3,131,117,115,8,123,5,105,2,124,106,104,7,125,4,107,1];
num_freedom_id = [(num-1).*3+1, (num-1).*3+2, (num-1).*3+3];
Stiff_total = Stiff;
Stiff_total(num_freedom_id,:) = [];
Stiff_total(:,num_freedom_id) = [];

Mass_total = Mass;
Mass_total(num_freedom_id,:) = [];
Mass_total(:,num_freedom_id) = [];

rank_Stiff_total = rank(Stiff_total);
rank_Mass_total = rank(Mass_total);

% s = load('C:\Users\win9\Documents\WeChat Files\wxid_1k6tvp1av82f22\FileStorage\File\2022-05\c3d20_matrix1.mat');
% ansys_stiff = s.Stiff;
% ansys_mass = s.Mass;

function [matrix,a_out] = GetMatrix(fid)
a=textscan(fid,'%d %d %d %d %f','Delimiter',',');%
a_max = max(a{5})
n = max(a{1})*3;
matrix = zeros(n,n);
num = 0;
for i = 1:length(a{1})
        
        matrix( 3*(a{1}(i)-1)+a{2}(i) , 3*(a{3}(i)-1)+a{4}(i) )=a{5}(i);
        matrix( 3*(a{3}(i)-1)+a{4}(i) , 3*(a{1}(i)-1)+a{2}(i) )=a{5}(i);
        if a{5}(i) ==1.0000e+36
            num = num+1;
            a_out(num,:)=[a{1}(i) a{2}(i) a{3}(i) a{4}(i)];
        end

end

end
