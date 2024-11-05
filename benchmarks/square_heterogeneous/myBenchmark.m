function mesh = myBenchmark()

global h edgTagToBC

edgTag = {1};
BC = {'ABC'};
edgTagToBC = containers.Map(edgTag,BC);

linkMsh = 'output/mesh.msh';
linkGeo = append(fileparts(mfilename('fullpath')),'/square_heterogeneous.geo');
system(['gmsh -2 ' linkGeo ' -o ' linkMsh ' -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh2D(linkMsh);

% -------------------------------------------------------------------------
% Physical coefficients
% -------------------------------------------------------------------------

% Physical parameters
global omega c1 c2 c3 rho1 rho2 rho3

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

    if y>0.25
        rho(tri) = rho1;
        c(tri) = c1;
    end
    if y>0 && y<0.25
        rho(tri) = rho2;
        c(tri) = c2;
    end
    if y>-0.25 && y<0
        if x>-0.25 && x<0.25
            rho(tri) = rho2;
            c(tri) = c2;
        else
            rho(tri) = rho3;
            c(tri) = c3;
        end
    end
    if y<-0.25
        rho(tri) = rho3;
        c(tri) = c3;
    end
    if y==0
        rho(tri) = rho2;
        c(tri) = c2;
    end
    if y==-0.25
        rho(tri) = rho3;
        c(tri) = c3;
    end

    eta(tri) = rho(tri) * c(tri);
    k(tri) = omega / c(tri);
end

end