function plotFieldExact2D(mesh)

global omega BCLeft BCRight

N = 2; % number of nodes per element

figure;
for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(tri,1));
    coord2 = mesh.coordV(mesh.listE(tri,2));
    x=linspace(coord1,coord2,N);
    [eta, ~, c] = physical_parameters_2D(mesh, tri);
    k = omega / c;
    [solP, ~] = mySol2D_heterogeneous(x,k,eta);
    plot(x,real(solP),'-b','LineWidth',1.5);
    hold on
    plot(x,imag(solP),'-r','LineWidth',1.5);
end
title('Exact solution');
legend('Real part', 'Imaginary part','FontSize',12);
grid on;
axis equal;

end