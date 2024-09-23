function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global k;

% Warning: This is not the reference solution. This is the incident wave
% field.

theta = 4*pi/10;
solU = -exp(1i*k*(cos(theta)*x+sin(theta)*y));
solF = 0*x;
solDx = -1i*k*cos(theta) * exp(1i*k*(cos(theta)*x+sin(theta)*y));
solDy = -1i*k*sin(theta) * exp(1i*k*(cos(theta)*x+sin(theta)*y));
solVx = 0*x;
solVy = 0*x;

end