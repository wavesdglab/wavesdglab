function [solU, solDx, solDy, solVx, solVy] = mySol(x,y)

global omega cAir cObj rhoAir rhoObj Rdisk

solU = scattDiskPenetrable_Solution(omega/cAir,omega/cObj,rhoAir,rhoObj,Rdisk,x,y,'scattered');
solDx = zeros(size(x));
solDy = zeros(size(x));
solVx = zeros(size(x));
solVy = zeros(size(x));

end