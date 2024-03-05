function mesh = myBenchmark()

global h edgTagToBC

% BCWest, BCNorth, BCEast, BCSouth
edgTag = {4};
BC = {'DIR'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'benchmarks/disk_heterogeneous/disk.msh';
linkGeo = 'benchmarks/disk_heterogeneous/disk.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end