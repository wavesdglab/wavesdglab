function mesh = myBenchmark()

global edgTagToBC h1 h2

% BCWest, BCNorth, BCEast, BCSouth
edgTag = {1, 2, 3, 4};
BC = {'ROB', 'ROB', 'ROB', 'ROB'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'benchmarks/open_heterogeneous/open_heterogeneous.msh';
linkGeo = 'benchmarks/open_heterogeneous/open_heterogeneous.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h1) ' -clmin ' num2str(h2)]);
mesh = readMesh2D(linkMsh);

% Physical parameters
global omega c1 c2 rho1 rho2

% Table for physical parameters
global rho c eta k
rho = zeros(mesh.numTri,1);
c = zeros(mesh.numTri,1);
eta = zeros(mesh.numTri,1);
k = zeros(mesh.numTri,1);

for tri = 1:mesh.numTri
    verTri = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(verTri(1),:);
    V2 = mesh.coord(verTri(2),:);
    V3 = mesh.coord(verTri(3),:);
    x = (V1(1,1)+V2(1,1)+V3(1,1))/3;

    if min(x) < 0.5
        rho(tri) = rho1;
        c(tri) = c1;
    else
        rho(tri) = rho2;
        c(tri) = c2;
    end

    eta(tri) = rho(tri) * c(tri);
    k(tri) = omega / c(tri);

end