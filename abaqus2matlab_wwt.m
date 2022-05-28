%%��abaqus�����ļ�����Ϊmatlab�����ļ�
%��ȡabaqus�Ľڵ�͵�Ԫ��Ϣ
%abaqus�ڵ���Ϣ ���� ��� xyz���� ��ά����
clc;clear;close all;
Dir = 'D:\SIMULIA\Temp';

hed = 'C3D20_cantilever_beam_matlab';                  %�ļ�ͷ
hed_Flag = 2;                                                            % 1-C3D8  2-C3D20
numnp = 321;                                                             %�ڵ����
numel = 40;                                                                 %��Ԫ����
ID(1:3,numnp)=0;
ID(1:3,[9,130,6,116,3,131,117,115,8,123,5,105,2,124,106,104,7,125,4,107,1])=1; 
%���Լ��

abaqus_inp_ID = ['\Job_',hed,'.inp'];%abaqus·��
matlab_in_ID = ['\Job_',hed,'.in'];%matlab·��

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

numnp = length(node_data{1});%�ڵ����
numel = length(element_data{1});%��Ԫ����

numeg = 1;%��Ԫ������
numem = 1;%��Ԫ����
%length(npar(:,1)) = numeg
if hed_Flag == 1
    npar(1,1) = 3;                                                              %��Ԫ���� 3-C3D8 4-C3D20
elseif hed_Flag == 2
    npar(1,1) = 4;
else
    disp('��Ԫ����npar1���');
end
npar(1,2) = numel;%����Ԫ���е�Ԫ����
npar(1,3) = 1;%����Ԫ���в�ͬ����/������������
    %length(material(:,1)) = npar(1,3)
    material(1,1) = 210000;                                              %E(1)   ����ģ��
    material(1,2) = 0.3;                                                     %MU(1)  ���ɱ�
    material(1,3) = 7.8e-9;                                                   %RHO(1) �ܶ� 1*10^3[kg/(mm^3)]
npar(1,4) = 2;                                                                  %��˹���ֽ״�

nlcase = 1;%�غɹ�����
    nload(1) = 16;%�����غɵĸ���
    %length(load_data(:,1)) = nload(1)
    load_data(1:nload(1),2) = [93 99 310 320 311 319 96 309 304 315 301 314 302 91 97 94];
    %�����غ����õĽڵ��
    load_data(1:nload(1),3) = 1.*ones(16,1);                                       %�غ����÷���
    load_data(1:nload(1),4) = [-833.33 -833.33 3333.3 3333.3 3333.3 3333.3 ...
        -1666.7 6666.7 -3333.3 -3333.3 -3333.3 -3333.3 -6666.7 833.33 833.33 1666.7];                                   %�غ�����ֵ
modex = 2;%���ģʽ


fileID_matlab = fopen([Dir,matlab_in_ID], 'w+');
%���������
fprintf(fileID_matlab,hed);
fprintf(fileID_matlab,'\n');
%���������
fprintf(fileID_matlab,' %4d %4d %4d %4d %4d\n', numnp,numeg,numem,nlcase,modex);
%����ڵ�����
for i = 1:numnp
    fprintf(fileID_matlab,' %4d %4d %4d %4d  %17.9E %17.9E %17.9E\n',i,ID(1,i),ID(2,i),ID(3,i),node_data{2}(i),node_data{3}(i),node_data{4}(i));
end
%����غ�����
for i = 1:nlcase
    fprintf(fileID_matlab,' %4d %4d\n',i,nload(i));
    for j = 1:nload(i)
        fprintf(fileID_matlab,' %4d %4d %9d\n',load_data(j,2),load_data(j,3),load_data(j,4));
    end
end
%��ά��Ԫ����
for i = 1:numeg
    %��Ԫ���������
    fprintf(fileID_matlab,' %4d %4d %4d %4d\n',npar(i,1),npar(i,2),npar(i,3),npar(i,4));
    %����/������������
    for j = 1:npar(i,3)
        fprintf(fileID_matlab,' %4d %9d %4.1f %9.1E\n',j,material(j,1),material(j,2),material(j,3));
    end
    %��Ԫ����
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