function mesh = myBenchmark()

<<<<<<< HEAD
global edgTagToBC h1 h2
=======
global h1 h2 edgTagToBC
>>>>>>> 51d25c5f96bc7ec1f21898009f162121a8e9812e

% BCWest, BCNorth, BCEast, BCSouth
edgTag = {1, 2, 3, 4};
BC = {'ROB', 'ROB', 'ROB', 'ROB'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'benchmarks/open_heterogeneous/open_heterogeneous.msh';
linkGeo = 'benchmarks/open_heterogeneous/open_heterogeneous.geo';
system(['gmsh -2 ' linkGeo ' -v 0 -o ' linkMsh ' -setnumber h1 ' num2str(h1) ' -setnumber h2 ' num2str(h2)]);
mesh = readMesh2D(linkMsh);

<<<<<<< HEAD
=======
% -------------------------------------------------------------------------
% Physical coefficients
% -------------------------------------------------------------------------

>>>>>>> 51d25c5f96bc7ec1f21898009f162121a8e9812e
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
<<<<<<< HEAD
=======
    y = (V1(1,2)+V2(1,2)+V3(1,2))/3;
>>>>>>> 51d25c5f96bc7ec1f21898009f162121a8e9812e

    if min(x) < 0.5
        rho(tri) = rho1;
        c(tri) = c1;
    else
        rho(tri) = rho2;
        c(tri) = c2;
    end

    eta(tri) = rho(tri) * c(tri);
    k(tri) = omega / c(tri);
<<<<<<< HEAD
=======
end
>>>>>>> 51d25c5f96bc7ec1f21898009f162121a8e9812e

end