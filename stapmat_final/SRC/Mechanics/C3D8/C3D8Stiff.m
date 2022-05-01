% % Called by GetStiff.m
% % 
% % - Call procedures:
% %        InitC3D8() % 定义每个单元上的节点数和每个节点的自由度
% %        ReadC3D8()
% % 

function C3D8Stiff()

% Init variables of the element
InitC3D8();

% Read Material and Elements
ReadC3D8();

fprintf('Solution phase ...\n\n');

% calculate addresses of diagonal elements
Addres();

end % end of Solid3DStiff

% ----------------------- Functions -----------------------------------

% Init parameters of C3D8 element
function InitC3D8()

global sdata;
sdata.NNODE = 8; % 一个单元上的节点数
sdata.NDOF = 3; % 每个节点的自由度

end
