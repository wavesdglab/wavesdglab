function [eta, rho, c] = physical_parameters_1D(mesh, e)

x_C = (mesh.coordV(1,e) +  mesh.coordV(1,e+1)) / 2;

if x_C < 1/2
    rho = 10;
    c = 1.1;
else
    rho = 0.1;
    c = 1.1;
end

eta = rho * c;

end