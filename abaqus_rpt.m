clc;clear;close all;

file_ID = 'D:\SIMULIA\Temp\abaqus-C3D8-fix4p-explicit-time20.rpt';
% file_ID = 'D:\SIMULIA\Temp\abaqus-C3D8-canbeam-explicit.rpt';
% result = importdata(file_ID,' ',19);
node_num = 8; %算例节点个数
node_plot = 6;

for i = 1:41
    fid= fopen(file_ID);
    result{i} =textscan(fid,'%d %f %f %f',node_num,'Delimiter',',','headerlines',19+(30+node_num )*(i-1));
    abaqus_rst(i) = result{i}{3}(node_plot);
    fclose(fid);
end

% node_u(:,[1 2 3])=result{i}(:,[2 3 4]);
% node_u(1,1:321)=str2num(data.textdata{14:334,2});
% for i = 1:length(node_u)
%     abaqus_rst_v1([(i-1)*3+1,(i-1)*3+2,(i-1)*3+3],1) = node_u(i,:);
% end

%去除固定约束自由度
bc_node= [9,6,3,8,5,2,7,4,1];
% bc_node= [1,3,5];
bc_length = length(bc_node);
% bc_node = [1,3,5,7];
% num_freedom_id = [(bc_node-1).*3+1, (bc_node-1).*3+2, (bc_node-1).*3+3];
% abaqus_rst = abaqus_rst_v1;
% abaqus_rst(num_freedom_id,:) = [];

% matlab_load = load('C:\Users\win9\Documents\WeChat Files\wxid_1k6tvp1av82f22\FileStorage\File\2022-05\C3D8_result-a_time0-1_step0.05.mat');
% matlab_rst = matlab_load.a((node_plot-1-bc_length)*3+1,:);

time = linspace(0,20,41);

figure;
plot(time,abaqus_rst,'LineWidth',2);
% hold on;
% plot(time,matlab_rst,'LineWidth',2);
% hold off;
set(gca,'FontSize',14);
% xlim([200 270]);
ylim([-1*10^(-3),0]);
% legend('abaqus','matlab','Location','southeast');
xlabel('time (s)');
ylabel('U8');
legend('abaqus');