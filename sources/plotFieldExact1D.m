function plotFieldExact1D(mesh)

global omega BCLeft BCRight

N = 2; % number of nodes per element

figure;
for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    x=linspace(coord1,coord2,N);
    [eta, ~, c] = physical_parameters_1D(mesh, e);
    k = omega / c;
    [solP, ~] = mySol1D_heterogeneous(x,k,eta);
    plot(x,real(solP),'-b');
    hold on
    plot(x,imag(solP),'-r');
end
title('Exact solution');
legend('Real part', 'Imaginary part','FontSize',15);
grid on;
axis equal;

end