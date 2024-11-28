function [solU, solDx, solDy, solVx, solVy] = mySol(x,y)

global k

[solU, solDx, solDy] = cavity(x,y,k);
solVx = solDx/(1i*k);
solVy = solDy/(1i*k);

end