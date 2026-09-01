function mesh = myBenchmark()

global h edgTagToBC

edgTag = {1, 2, 3, 4, 5, 6, 7};
BC = {'ABC', 'ABC', 'NEU0', 'NEU0', 'ABC', 'ABC', 'ABC'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'output/mesh.msh';

% Mesh with interface
linkGeo = 'benchmarks/aeroacou_hydrointerf/mesh.geo';

% Mesh without interface
% linkGeo = 'benchmarks/aeroacou_hydrointerf/mesh2.geo';

system(['gmsh -2 ' linkGeo ' -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end