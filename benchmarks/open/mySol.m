function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global k theta

solU  = exp(1i*k*(cos(theta)*x+sin(theta)*y));
solF  = 0*x;
solDx = 1i*k*cos(theta) * solU;
solDy = 1i*k*sin(theta) * solU;
solVx = cos(theta) * solU;
solVy = sin(theta) * solU;

end
