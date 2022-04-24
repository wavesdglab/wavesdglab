close all;
clear all;

headers1D()

x = -1:0.01:1;
degree = 10;

% Shape functions

figure(1);
val = functionsShapeLIN(x,degree);
hold off
plot(x,val,'Linewidth',2);
ylim([-1 2]);
legend('$\phi_0$','$\phi_1$','$\phi_2$','$\phi_3$','$\phi_4$','Interpreter','latex','FontSize',16);
xlabel('$u$','Interpreter','latex','FontSize',16);
title('Shape Functions');
box on

Q = 10;
[nodes, weights] = quadratureGaussLIN(Q);
shapeFunc = functionsShape1D(nodes,degree);
for i=1:degree+1
    L2norm(i) = weights(:)' * abs(shapeFunc(:,i));
end
L2norm = sqrt(L2norm);
L2norm'

% Derivative of the shape functions

figure(2);
val = functionsShapeDerLIN(x,degree);
plot(x,val,'Linewidth',2);
ylim([-2 2]);
legend('$d\phi_0$','$d\phi_1$','$d\phi_2$','$d\phi_3$','$d\phi_4$','Interpreter','latex','FontSize',16);
xlabel('$u$','Interpreter','latex','FontSize',16);
title('Shape Functions (derivative)');
box on
% 

Q = 10;
[nodes, weights] = quadratureGaussLIN(Q);
shapeFunc = functionsShape1D(nodes,degree);
for i=1:degree+1
    L2norm(i) = weights(:)' * abs(shapeFunc(:,i));
end
L2norm = sqrt(L2norm);
L2norm'
