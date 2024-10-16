function mesh = myBenchmark()

global edgTagToBC h1 h2

% BCWest, BCNorth, BCEast, BCSouth
edgTag = {4};
BC = {'DIR'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'benchmarks/disk_heterogeneous/disk.msh';
linkGeo = 'benchmarks/disk_heterogeneous/disk.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -clmax ' num2str(h1) ' -clmin ' num2str(h2)]);
mesh = readMesh2D(linkMsh);

% -------------------------------------------------------------------------
% Physical coefficients
% -------------------------------------------------------------------------

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
    y = (V1(1,2)+V2(1,2)+V3(1,2))/3;

    R1 = 0.25;
    if x^2 + y^2 < R1^2
        rho(tri) = rho1;
        c(tri) = c1;
    else
        rho(tri) = rho2;
        c(tri) = c2;
    end

    eta(tri) = rho(tri) * c(tri);
    k(tri) = omega / c(tri);
end

end