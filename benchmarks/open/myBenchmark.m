function mesh = myBenchmark()

global h edgTagToBC theta

theta = pi/4;

% BCWest, BCNorth, BCEast, BCSouth
edgTag = {1, 2, 3, 4};
BC = {'ROB', 'ROB', 'ROB', 'ROB'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'output/mesh.msh';
linkGeo = append(fileparts(mfilename('fullpath')),'/open.geo');
system(['gmsh -2 ' linkGeo ' -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end