function [eta, rho, c] = physical_parameters(mesh, tri)

verTri = mesh.mapTriToVer(tri,:);
V1 = mesh.coord(verTri(1),:);
V2 = mesh.coord(verTri(2),:);
V3 = mesh.coord(verTri(3),:);
x_C = (V1(1,1)+V2(1,1)+V3(1,1))/3;
% y_C = (V1(1,2)+V2(1,2)+V3(1,2))/3;

% % Cavity / Open
% if (x_C<0.5)
%     rho = 1;
% %     eta = 1;
%     c = 1;
% else
%     rho = 1;
% %     eta = 8;
%     c = 1.5;
% end
% eta = rho*c;
% % rho = eta / c;

% Waveguide
if (x_C<2)
    rho = 1;
%     eta = 1;
    c = 0.5;
else
    rho = 1;
%     eta = 8;
    c = 1;
end
eta = rho*c;
% rho = eta / c;

end