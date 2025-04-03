function mesh = myBenchmark()

global h edgTagToBC M hmin
global pntSouTag pntSouVal

edgTag = {3};
BC = {'ABC'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'output/mesh.msh';
linkGeo = 'benchmarks/disk_convected/disk.geo';

hmin = h/4;

% system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -setnumber h ' num2str(h) ' -setnumber hmin ' num2str(hmin) ' -setnumber M ' num2str(M)]);

system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ' -setnumber M ' num2str(M)]);
% system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

pntSouTag = 1;
pntSouVal = 1;
vertSou = mesh.mapPntToVer(mesh.tagPntFile == pntSouTag);
xSou = mesh.coord(vertSou,1);
ySou = mesh.coord(vertSou,2);

end