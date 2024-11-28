function [solU, solDx, solDy, solVx, solVy] = mySol(x,y)

global k Rdisk

solU = solScattPlaneWaveHard(k,Rdisk,x,y);
solDx = -1i*k*exp(1i*k*x);
solDy = 0*x;
solVx = 0*x;
solVy = 0*x;

end