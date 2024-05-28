function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global k Rdisk

solU = solScattPlaneWaveHard(k,Rdisk,x,y);
solF = 0*x;
solDx = -1i*k*exp(1i*k*x);
solDy = 0*x;
solVx = 0*x;
solVy = 0*x;

end