%%对比abaqus和matlab计算结果
%1.节点位移 2.节点应力 3.节点应变
clc;clear;close all;
Dir = 'D:\ProgramData\Dynamic_Mechanic\Computational-dynamics-2022\stapmat_final\Data\';
hed = 'C3D20_ne1_P_wwt';                              %文件头
abaqus_out_ID = [hed,'_abaqus','.rpt'];            %abaqus路径
matlab_out_ID_1 = [hed,'.OUT'];                      %matlab路径
matlab_out_ID_2 = [hed,'_curv','.DAT'];            %matlab路径
