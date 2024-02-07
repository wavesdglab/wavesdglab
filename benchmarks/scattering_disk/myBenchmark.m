function mesh = myBenchmark()

global h edgTagToBC
global LdomX LdomY LpmlX LpmlY Rdisk

if(isempty(LdomX)) LdomX = 1.2; end
if(isempty(LdomY)) LdomY = 1.2; end
if(isempty(LpmlX)) LpmlX = 0.3; end
if(isempty(LpmlY)) LpmlY = 0.3; end
if(isempty(Rdisk)) Rdisk = 1; end

% BCPML, BCObstacle
edgTag = {1, 2};
BC = {'DIR0', 'NEU'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'benchmarks/scattering_disk/scatteringPML.msh';
linkGeo = 'benchmarks/scattering_disk/scatteringPML.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ' -setnumber LdomX ' num2str(LdomX) ' -setnumber LdomY ' num2str(LdomY) ' -setnumber LpmlX ' num2str(LpmlX) ' -setnumber LpmlY ' num2str(LpmlY) ' -setnumber Rdisk ' num2str(Rdisk)]);
mesh = readMesh2D(linkMsh);

end