function [solU, solDx, solDy, solVx, solVy] = mySol(x,y)

global omega cAir cObj rhoAir rhoObj Rdisk

solU = scattDiskPenetrable_Solution(omega/cAir,omega/cObj,rhoAir,rhoObj,Rdisk,x,y);
solDx = 0;
solDy = 0;
solVx = 0;
solVy = 0;

end