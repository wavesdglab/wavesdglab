function mesh = myBenchmark()

global h edgTagToBC
global L L_PML

% BCPML, BCObstacle
edgTag = {1, 2};
BC = {'DIR', 'NEU'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'benchmarks/scattering_square/square_cavity.msh';
linkGeo = 'benchmarks/scattering_square/square_cavity.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ' -setnumber L ' num2str(L) ' -setnumber L_PML ' num2str(L_PML)]);
mesh = readMesh2D(linkMsh);

end