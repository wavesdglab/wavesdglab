function mesh = myBenchmark()

global h edgTagToBC

% BCWest, BCNorth, BCEast, BCSouth
edgTag = {3};
BC = {'ABC'};
% BC = {'ROB'};
edgTagToBC = containers.Map(edgTag,BC);

% linkMsh = 'benchmarks/waveguide_convected/waveguide.msh';
linkMsh = 'output/mesh.msh';
linkGeo = 'benchmarks/disk_convected/disk.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end