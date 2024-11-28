function [solU, solDx, solDy, solVx, solVy] = mySol(x,y)

global k theta

solU  = exp(1i*k*(cos(theta)*x+sin(theta)*y));
solDx = 1i*k*cos(theta) * solU;
solDy = 1i*k*sin(theta) * solU;
solVx = cos(theta) * solU;
solVy = sin(theta) * solU;

end
