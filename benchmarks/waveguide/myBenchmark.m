function mesh = myBenchmark()

global h edgTagToBC

% BCWest, BCNorth, BCEast, BCSouth
edgTag = {1, 2, 3, 4};
BC = {'DIR0', 'DIR0', 'ROB', 'DIR0'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'output/mesh.msh';
linkGeo = append(fileparts(mfilename('fullpath')),'/waveguide.geo');
system(['gmsh -2 ' linkGeo ' -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end