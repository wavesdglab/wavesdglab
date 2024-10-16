function mesh = myBenchmark()

global h edgTagToBC
global LdomX LdomY LpmlX LpmlY

if(isempty(LdomX)) LdomX = 1.1; end
if(isempty(LdomY)) LdomY = 1.1; end
if(isempty(LpmlX)) LpmlX = 0.1; end
if(isempty(LpmlY)) LpmlY = 0.1; end

% BCPML, BCObstacle
edgTag = {1, 2};
BC = {'DIR0', 'DIR'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'output/mesh.msh';
linkGeo = 'benchmarks/scattering_openCavity_DIR/scattering_openCavity.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h) ' -setnumber LdomX ' num2str(LdomX) ' -setnumber LdomY ' num2str(LdomY) ' -setnumber LpmlX ' num2str(LpmlX) ' -setnumber LpmlY ' num2str(LpmlY)]);
mesh = readMesh2D(linkMsh);

end