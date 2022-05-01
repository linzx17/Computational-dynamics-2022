function vtkwrite(name,i)
global cdata;
global sdata;
ID = sdata.ID; 
NUME = cdata.NPAR(2);NUMNP = cdata.NUMNP;%输出使用的数据
DIS = sdata.DIS;
VV = sdata.V;
STRAIN = sdata.STRAIN;
STRESS = sdata.STRESS;
D = zeros(NUMNP, 3, 'double');
V = zeros(NUMNP, 3, 'double');
filename = strcat('.\Data\',cdata.HED);
midname = strcat(filename,'\mid.vtk');
filename = strcat(filename,'\',name,'N=',num2str(i),'.vtk');
copyfile (midname, filename)
fid = fopen(filename, 'a');
fprintf(fid,'POINT_DATA %d\n',NUMNP);
for II = 1:NUMNP
    if (ID(1, II) ~= 0)
        D(II,1) = DIS(ID(1, II));
        V(II,1) = VV(ID(1, II));
    end
    if (ID(2, II) ~= 0)
        D(II,2) = DIS(ID(2, II));
        V(II,2) = VV(ID(2, II));
    end
end
fprintf(fid,'VECTORS dis double\n');
dlmwrite(filename,D,'-append',...
    'delimiter',' ')
fprintf(fid,'VECTORS vec double\n');
dlmwrite(filename,V,'-append',...
    'delimiter',' ')


fprintf(fid,'CELL_DATA %d\n',NUME);
fprintf(fid,'SCALARS SRR double 1\n');
fprintf(fid,'LOOKUP_TABLE srr\n');
fprintf(fid,'%e\n',STRESS(:,1));

fprintf(fid,'SCALARS SZZ double 1\n');
fprintf(fid,'LOOKUP_TABLE szz\n');
fprintf(fid,'%e\n',STRESS(:,2));

fprintf(fid,'SCALARS SRZ double 1\n');
fprintf(fid,'LOOKUP_TABLE srz\n');
fprintf(fid,'%e\n',STRESS(:,3));

fprintf(fid,'SCALARS STHETA double 1\n');
fprintf(fid,'LOOKUP_TABLE stheta\n');
fprintf(fid,'%e\n',STRESS(:,4));

fprintf(fid,'SCALARS Mises double 1\n');
fprintf(fid,'LOOKUP_TABLE stheta\n');
fprintf(fid,'%e\n',STRESS(:,5));

fprintf(fid,'SCALARS DxRR double 1\n');
fprintf(fid,'LOOKUP_TABLE dxrr\n');
fprintf(fid,'%e\n',STRAIN(:,1));

fprintf(fid,'SCALARS DxZZ double 1\n');
fprintf(fid,'LOOKUP_TABLE dxzz\n');
fprintf(fid,'%e\n',STRAIN(:,2));

fprintf(fid,'SCALARS DxRZ double 1\n');
fprintf(fid,'LOOKUP_TABLE dxrz\n');
fprintf(fid,'%e\n',STRAIN(:,3));

fprintf(fid,'SCALARS DxTHETA double 1\n');
fprintf(fid,'LOOKUP_TABLE dxtheta\n');
fprintf(fid,'%e\n',STRAIN(:,4));
end