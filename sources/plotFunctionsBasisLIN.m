% This script plots the basis functions for 5 different bases

function plotFunctionsBasisLIN()

x = -1:0.01:1;
nodes = [-1 -0.5 0 0.5 1];
degree = 4;
N = degree+1;

% Monomials

figure(1);
hold off
plot(x,x.^0,'Linewidth',2);
hold on
plot(x,x.^1,'Linewidth',2);
plot(x,x.^2,'Linewidth',2);
plot(x,x.^3,'Linewidth',2);
plot(x,x.^4,'Linewidth',2);
legend('$m_0$','$m_1$','$m_2$','$m_3$','$m_4$','Interpreter','latex','FontSize',16);
xlabel('$u$','Interpreter','latex','FontSize',16);
title('Monomials');
box on

% Lagrange functions

figure(2);
val = functionsLagrange(x,nodes);
plot(x,val,'Linewidth',2);
ylim([-1 2]);
legend('$\ell_0$','$\ell_1$','$\ell_2$','$\ell_3$','$\ell_4$','Interpreter','latex','FontSize',16);
xlabel('$u$','Interpreter','latex','FontSize',16);
title('Lagrange functions');
box on

% Legendre functions

figure(3);
val = functionsLegendre(x,degree);
plot(x,val,'Linewidth',2);
%ylim([-1 1]);
legend('$L_0$','$L_1$','$L_2$','$L_3$','$L_4$','Interpreter','latex','FontSize',16);
xlabel('$u$','Interpreter','latex','FontSize',16);
title('Legendre functions');
box on

% Lobbato functions

figure(4);
val = functionsLobbato(x,degree);
plot(x,val,'Linewidth',2);
ylim([-1 1]);
legend('$\phi_0$','$\phi_1$','$\phi_2$','$\phi_3$','$\phi_4$','Interpreter','latex','FontSize',16);
xlabel('$u$','Interpreter','latex','FontSize',16);
title('Lobbato functions');
box on

% Bernstein functions

figure(5);
val = functionsBernstein(x,degree);
plot(x,val,'Linewidth',2);
ylim([-1 1]);
legend('$B_0$','$B_1$','$B_2$','$B_3$','$B_4$','Interpreter','latex','FontSize',16);
xlabel('$u$','Interpreter','latex','FontSize',16);
title('Bernstein functions');
box on

end