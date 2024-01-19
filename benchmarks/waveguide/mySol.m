function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global k;

L = 1.;
L = 4.;
theta = 30.*(pi/180.);
[solU, solDx, solDy] = waveguide(x,y,k,L,theta);
solF = 0*x;
solVx = solDx/(1i*k);
solVy = solDy/(1i*k);

end