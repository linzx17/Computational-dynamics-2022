%%��abaqus�����ļ�����Ϊmatlab�����ļ�
%��ȡabaqus�Ľڵ�͵�Ԫ��Ϣ
%abaqus�ڵ���Ϣ ���� ��� xyz���� ��ά����
clc;clear;close all;
Dir = 'D:\SIMULIA\Temp';

hed = 'C3D20-cantilever-beam-element320';                  %�ļ�ͷ
hed_Flag = 2;                                                            % 1-C3D8  2-C3D20
numnp = 1865;                                                             %�ڵ����
numel = 320;                                                                 %��Ԫ����
numeg = 1;                                                                  %��Ԫ������
numem = 320;                                                                 %��Ԫ����
ID(1:3,numnp)=0;
bc_node = [25 628 20 605 15 582 10 558 5 ...
                    629 606 583 559 557 ...
                    24 623 19 600 14 577 9 550 4 ...
                    624 601 578 551 549 ...
                    23 618 18 595 13 572 8 542 3 ...
                    619 596 573 543 541 ...
                    22 611 17 588 12 565 7 531 2 ...
                    612 589 566 532 530 ...
                    21 613 16 590 11 567 6 533 1];
ID(1:3,bc_node)=1; 
%���Լ��

abaqus_inp_ID = ['\Job-',hed,'-matlab.inp'];%abaqus·��
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
npar(1,4) = 3;                                                                  %��˹���ֽ״�

M = 100000;%[N��mm]
L = 5;%[mm]
S = 5*10;%[mm^2]
q_surfload = M/S/L;%[N/mm^2]
a = 2.5;%[mm]
q_node = a^2*q_surfload*(-1/12);
q_egde = a^2*q_surfload*(1/3);
q_load = zeros(65,1);%�кŽڵ��ţ�ֵ �غɴ�С
num=[ ...
    1 3 15 17 2 10 11 16 ;
    3 5 17 19 4 11 12 18 ;
    5 7 19 21 6 12 13 20 ;
    7 9 21 23 8 13 14 22 ;
    15 17 29 31 16 24 25 30 ;
    17 19 31 33 18  25 26 32 ;
    19 21 33 35 20 26  27 34 ;
    21 23 35 37 22 27 28 36 ;
    29 31 43 45 30 38 39 44 ;
    31 33 45 47 32 39 40 46 ;
    33 35 47 49 34 40 41 48 ;
    35 37 49 51 36 41 42 50 ;
    43 45 57 59 44 52 53 58 ;
    45 47 59 61 46 53 54 60 ;
    47 49 61 63 48 54 55 62 ;
    49 51 63 65 50 55 56 64 ];
for i = 1:16
    if i <9
        q_load(num(i,1:4),1) = q_load(num(i,1:4),1) + q_node.*ones(4,1);
        q_load(num(i,5:8),1) = q_load(num(i,5:8),1) + q_egde.*ones(4,1);
    else
        q_load(num(i,1:4),1) = q_load(num(i,1:4),1) - q_node.*ones(4,1);
        q_load(num(i,5:8),1) = q_load(num(i,5:8),1) - q_egde.*ones(4,1);
    end
end
load_node = [505 1820 510 1836 515 1850 520 1864 525 ...
                      1821 1819 1835 1849 1863 ...
                      504 1815 509 1833 514 1847 519 1861 524 ...
                      1816 1814 1832 1846 1860 ...
                      503 1810 508 1830 513 1844 518 1858 523 ...
                      1811 1806 1829 1643 1857 ...
                      502 1803 507 1826 512 1840 517 1854 522 ...
                      1804 1802 1825 1839 1853 ...
                      501 1801 506 1824 511 1838 516 1852 521];

nlcase = 1;%�غɹ�����
    nload(1) = 65;%�����غɵĸ���
    %length(load_data(:,1)) = nload(1)
    load_data(1:nload(1),2) = load_node;
    %�����غ����õĽڵ��
    load_data(1:nload(1),3) = 1.*ones(nload(1),1);                                       %�غ����÷���
    load_data(1:nload(1),4) =  q_load;                              %�غ�����ֵ
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
                ,element_data{2}(k),element_data{3}(k),element_data{4}(k),element_data{5}(k)...
                ,element_data{6}(k),element_data{7}(k),element_data{8}(k),element_data{9}(k)...
                ,numeg...
                );
        elseif length(element_data)-1 == 20 %C3D20
            fprintf(fileID_matlab,' %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n',k...
                ,element_data{2}(k),element_data{3}(k),element_data{4}(k),element_data{5}(k)...
                ,element_data{6}(k),element_data{7}(k),element_data{8}(k),element_data{9}(k)...
                ,element_data{10}(k),element_data{11}(k),element_data{12}(k),element_data{13}(k)...
                ,element_data{14}(k),element_data{15}(k),element_data{16}(k),element_data{17}(k)...
                ,element_data{18}(k),element_data{19}(k),element_data{20}(k),element_data{21}(k)...
                ,numeg...
                );
        else
            fprintf(fileID_matlab,'Unknown solid element type!\n');
            fprintf(fileID_matlab,'Check element_data!\n');
        end
    end
end
fprintf(fileID_matlab,'stop'); 