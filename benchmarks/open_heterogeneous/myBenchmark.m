function mesh = myBenchmark()

global h edgTagToBC

% BCWest, BCNorth, BCEast, BCSouth
edgTag = {1, 2, 3, 4};
BC = {'ROB', 'ROB', 'ROB', 'ROB'};
edgTagToBC = containers.Map(edgTag,BC);

% linkMsh = 'benchmarks/open_heterogeneous/open_heterogeneous.msh';
% linkGeo = 'benchmarks/open_heterogeneous/open_heterogeneous.geo';
% system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
% mesh = readMesh2D(linkMsh);

h1=1;
h2=0;

linkMsh = 'benchmarks/open_heterogeneous/open_heterogeneous.msh';
linkGeo = 'benchmarks/open_heterogeneous/open_heterogeneous.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h1) ' -clmin ' num2str(h2)]);
mesh = readMesh2D(linkMsh);

end