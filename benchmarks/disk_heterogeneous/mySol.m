function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global eta1 eta2 k1 k2 c1 c2

R1 = 0.25;
R2 = 0.5;
f = 1;

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

r = sqrt(x.^2+y.^2);

if (max(r)<R1)
    solU = A_1 * besselj(0,k1*r) - f/k1^2;
    solDx = - k1 * A_1 * besselj(1,k1*r) .* x ./ r;
    solDy = - k1 * A_1 * besselj(1,k1*r) .* y ./ r;
    solF = f + 0*x + 0*y;
    solVx = 1i/eta1 * A_1 * besselj(1,k1*r) .* x ./ r;
    solVy = 1i/eta1 * A_1 * besselj(1,k1*r) .* y ./ r;
else
    solU = A_2 * besselj(0,k2*r) + B_2 * bessely(0,k2*r) - f/k2^2;
    solDx = - k2 * (A_2 * besselj(1,k2*r) + B_2 * bessely(1,k2*r)) .* x ./ r;
    solDy = - k2 * (A_2 * besselj(1,k2*r) + B_2 * bessely(1,k2*r)) .* y ./ r;
    solF = f + 0*x + 0*y;
    solVx = 1i/eta2 * (A_2 * besselj(1,k2*r) + B_2 * bessely(1,k2*r)) .* x ./ r;
    solVy = 1i/eta2 * (A_2 * besselj(1,k2*r) + B_2 * bessely(1,k2*r)) .* y ./ r;
end

%% Homogeneous medium

% R1 = 0.25;
% R2 = 0.5;
% f = 1;
% 
% A = f/(k1^2*besselj(0,k1*R2));
% 
% r = sqrt(x.^2+y.^2);
% 
% solU = A * besselj(0,k1*r) - f/(k1^2);            
% solDx = - k1 * A * besselj(1,k1*r) .* x ./ r;      
% solDy = - k1 * A * besselj(1,k1*r) .* y ./ r;       
% solF = f + 0*x + 0*y;                               
% solVx = 1i/eta1 * A * besselj(1,k1*r) .* x ./ r;    
% solVy = 1i/eta1 * A * besselj(1,k1*r) .* y ./ r;    

end