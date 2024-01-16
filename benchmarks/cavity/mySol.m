function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global k

[solU, solDx, solDy] = cavity(x,y,k);
solF = 0*x+1;
solVx = solDx/(1i*k);
solVy = solDy/(1i*k);

end