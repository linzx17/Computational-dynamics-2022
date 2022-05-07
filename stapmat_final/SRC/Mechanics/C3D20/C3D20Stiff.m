% % Called by GetStiff.m
% % 
% % - Call procedures:
% %        InitC3D20() % 定义每个单元上的节点数和每个节点的自由度
% %        ReadC3D20()
% % 

function C3D20Stiff()

% Init variables of the element
InitC3D20();

% Read Material and Elements
ReadC3D8();

fprintf('Solution phase ...\n\n');

% calculate addresses of diagonal elements
Addres();

% Data check Or Solve
global cdata;
if (cdata.MODEX == 0) 
    cdata.TIM(3,:) = clock;
    cdata.TIM(4,:) = clock;
    cdata.TIM(5,:) = clock;
    return; 
end

% Assemble structure stiffness matrix
Assemble();

end % end of C3D20Stiff()

% ----------------------- Functions -----------------------------------

% Init parameters of C3D8 element
function InitC3D20()

global sdata;
sdata.NNODE = 20; % 一个单元上的节点数
sdata.NDOF = 3; % 每个节点的自由度

end % end of function InitC3D20()


% Assemble structure stiffness matrix
function Assemble()

global sdata;
global cdata;
sdata.STIFF = zeros(sdata.NWK, 1, 'double');

NUME = sdata.NUME; MATP = sdata.MATP; XYZ = sdata.XYZ; 
E = sdata.E; nu = sdata.nu; NGaussian = sdata.NGaussian; LM = sdata.LM;


