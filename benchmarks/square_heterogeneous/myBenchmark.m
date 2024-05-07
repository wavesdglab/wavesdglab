function mesh = myBenchmark()

global h edgTagToBC

edgTag = {1};
BC = {'ABC'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'benchmarks/square_heterogeneous/square_heterogeneous.msh';
linkGeo = 'benchmarks/square_heterogeneous/square_heterogeneous.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end