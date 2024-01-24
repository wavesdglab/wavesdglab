function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global eta1 k1

L = 4.;
theta = 45.*(pi/180.);
[solU, solDx, solDy] = waveguide(x,y,k1,L,theta);
solF = 0*x;
solVx = solDx/(1i*eta1*k1);
solVy = solDy/(1i*eta1*k1);

end