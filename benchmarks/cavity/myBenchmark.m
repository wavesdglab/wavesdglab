function mesh = myBenchmark()

global h edgTagToBC

% BCWest, BCNorth, BCEast, BCSouth
edgTag = {1, 2, 3, 4};
BC = {'DIR0', 'DIR0', 'DIR0', 'DIR0'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'benchmarks/cavity/cavity.msh';
linkGeo = 'benchmarks/cavity/cavity.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end