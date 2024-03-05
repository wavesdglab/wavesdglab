function setParameters(mesh)

% Physical parameters
global omega eta1 eta2 k1 k2 c1 c2 rho1 rho2

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
    if min(x) < 0.5
        rho(tri) = rho1;
        c(tri) = c1;
    else
        rho(tri) = rho2;
        c(tri) = c2;
    end
%     if x^2 + y^2 < 0.25
%         rho(tri) = rho1;
%         c(tri) = c1;
%     else
%         rho(tri) = rho2;
%         c(tri) = c2;
%     end

    eta(tri) = rho(tri) * c(tri);
    k(tri) = omega / c(tri);
end