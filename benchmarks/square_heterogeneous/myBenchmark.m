function mesh = myBenchmark()

global h edgTagToBC

% BCWest, BCNorth, BCEast, BCSouth
edgTag = {1, 2, 3, 4};
BC = {'ABC', 'ABC', 'ABC', 'ABC'};
edgTagToBC = containers.Map(edgTag,BC)

linkMsh = 'benchmarks/square_heterogeneous/square_heterogeneous.msh';
linkGeo = 'benchmarks/square_heterogeneous/square_heterogeneous.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end