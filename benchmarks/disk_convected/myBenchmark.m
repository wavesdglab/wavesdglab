function mesh = myBenchmark()

global h edgTagToBC xSou ySou c v0d
global pntSouTag pntSouVal

edgTag = {3};
BC = {'ABC'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'output/mesh.msh';
linkGeo = 'benchmarks/disk_convected/disk.geo';

xSou = -v0d/c;
ySou = 0;
system(['gmsh -2 ' linkGeo ' -o ' linkMsh ' -setnumber xSou ' num2str(xSou) ' -setnumber ySou ' num2str(ySou) ' -setnumber h ' num2str(h)]);
mesh = readMesh2D(linkMsh);

pntSouTag = 1;
pntSouVal = 1;
%vertSou = mesh.mapPntToVer(mesh.tagPntFile == pntSouTag);
%xSou = mesh.coord(vertSou,1);
%ySou = mesh.coord(vertSou,2);

end