function [eta, rho, c] = physical_parameters(mesh, tri)

verTri = mesh.mapTriToVer(tri,:);
V1 = mesh.coord(verTri(1),:);
V2 = mesh.coord(verTri(2),:);
V3 = mesh.coord(verTri(3),:);
x_C = (V1(1,1)+V2(1,1)+V3(1,1))/3;
% y_C = (V1(1,2)+V2(1,2)+V3(1,2))/3;

% Cavity
if (x_C>0 && x_C<0.5)
    rho = 0.9;
    c = 0.8;
end
if (x_C>0.5 && x_C<1)
    rho = 0.9;
    c = 0.9;
end
eta = rho*c;

end