function mesh = myBenchmark()

global h edgTagToBC

% BCWest, BCNorth, BCEast, BCSouth
edgTag = {1, 2, 3, 4};
% BC = {'ROB', 'NEU0', 'ABC', 'NEU0'};
% BC = {'DIR', 'NEU0', 'ABC', 'NEU0'};
% BC = {'DIR', 'NEU0', 'DIR', 'NEU0'};
% BC = {'DIR', 'DIR', 'DIR', 'DIR'};

BC = {'ROB', 'ROB', 'ROB', 'ROB'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'benchmarks/waveguide_convected/waveguide.msh';
linkGeo = 'benchmarks/waveguide_convected/waveguide.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

end