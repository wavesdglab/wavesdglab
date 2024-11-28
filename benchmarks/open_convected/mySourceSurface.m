function [souU, souDx, souDy, souVx, souVy] = mySourceSurface(x,y)

global M rho c omega theta;

souU = zeros(size(x));
souDx = zeros(size(x));
souDy = zeros(size(x));
souVx = zeros(size(x));
souVy = zeros(size(x));

souU = (1+M)/2*exp(1i*omega/(c*(1+M)).*(x*cos(theta)+y*sin(theta)));

souVx = (1+M)/(2*rho*c)*exp(1i*omega/(c*(1+M)).*(x*cos(theta)+y*sin(theta)))*cos(theta);
souVy = (1+M)/(2*rho*c)*exp(1i*omega/(c*(1+M)).*(x*cos(theta)+y*sin(theta)))*sin(theta);

end