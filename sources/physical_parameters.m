function [eta, rho, c] = physical_parameters(mesh, tri)

verTri = mesh.mapTriToVer(tri,:);
V1 = mesh.coord(verTri(1),:);
V2 = mesh.coord(verTri(2),:);
V3 = mesh.coord(verTri(3),:);
x_C = (V1(1,1)+V2(1,1)+V3(1,1))/3;
% y_C = (V1(1,2)+V2(1,2)+V3(1,2))/3;

% Cavity / Open
if (x_C<0.5)
%     rho = 100;
    eta = 1;
    c = 1;
else
%     rho = 100;
    eta = 1;
    c = 0.8;
end
% eta = rho*c;
rho = eta / c;

end