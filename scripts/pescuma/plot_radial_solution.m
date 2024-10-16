function plot_radial_solution(k1,k2,eta1,eta2,R1,R2,f)

mat_11 = 0;
mat_12 = besselj(0,k2*R2); 
mat_13 = bessely(0,k2*R2);

mat_21 = besselj(0,k1*R1);
mat_22 = -besselj(0,k2*R1);
mat_23 = -bessely(0,k2*R1);

mat_31 = besselj(1,k1*R1)/eta1;
mat_32 = -besselj(1,k2*R1)/eta2;
mat_33 = -bessely(1,k2*R1)/eta2;

rhs_1 = f/k2^2;
rhs_2 = f/k1^2 - f/k2^2;
rhs_3 = 0;

mat = [mat_11 mat_12 mat_13;
       mat_21 mat_22 mat_23;
       mat_31 mat_32 mat_33];
rhs = [rhs_1;
       rhs_2;
       rhs_3];

coeff = mat \ rhs;

A_1 = coeff(1,1);
A_2 = coeff(2,1);
B_2 = coeff(3,1);

r1 = 0:0.0001:R1;
r2 = R1:0.0001:R2;

solU1 = A_1 * besselj(0,k1*r1) - f/k1^2;
solU2 = A_2 * besselj(0,k2*r2) + B_2 * bessely(0,k2*r2) - f/k2^2;

solU = [solU1 solU2];
r = [r1 r2];

figure
plot(r, solU, '-r');
grid on;
legend('p(r)')

% rezu1A = ["r" "solU"];
% rezu2A = [r', solU'];
% name = sprintf('output/plot_radial_sol_k1%g_k2%g.csv', k1, k2);
% writematrix([rezu1A ; rezu2A], name, 'Delimiter', 'semi');

det(mat)

end
