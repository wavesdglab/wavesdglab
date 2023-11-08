function [eta, rho, c] = physical_parameters_1D(mesh, e)

global omega

x_C = (mesh.coordV(1,e) +  mesh.coordV(1,e+1)) / 2;

% if x_C < 1/2
%     k = 0.01;
%     eta = 1;
% else
%     k = 0.5;
%     eta = 10;
% end
% 
% c = omega / k;
% rho = eta / c;

if x_C < 1/2
    rho = 0.1; 
    c = 2;
else
    rho = 0.1;
    c = 0.4;
end

eta = rho * c;

end