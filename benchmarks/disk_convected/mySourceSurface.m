function [souU, souDx, souDy, souVx, souVy] = mySourceSurface(x,y)

global M rho c omega theta phi;

xS=-M;
yS=0;

Theta = [cos(theta), sin(theta)];
Phi   = [cos(phi), sin(phi)];
souU  = zeros(size(x));
souDx = zeros(size(x));
souDy = zeros(size(x));
souVx = zeros(size(x));
souVy = zeros(size(x));

% souU = exp(1i*omega/(c*(1+M*(Theta*Phi'))).*(x*cos(phi)+y*sin(phi)));
% souVx = cos(phi)/(rho*c)*exp(1i*omega/(c*(1+M*(Theta*Phi'))).*(x*cos(phi)+y*sin(phi)));
% souVy = sin(phi)/(rho*c)*exp(1i*omega/(c*(1+M*(Theta*Phi'))).*(x*cos(phi)+y*sin(phi)));

% souU = 1i / (4*c^2*sqrt(1-M^2)) * besselh(0,omega/c*sqrt((x-xS).^2+(1-M^2)*(y-yS).^2)/(1-M^2)) .* exp(-1i*M/(1-M^2)*omega/c.*(x-xS));
souU = omega / (4*c^2*sqrt(1-M^2)) * besselh(0,omega/c*sqrt((x-xS).^2+(1-M^2)*(y-yS).^2)/(1-M^2)) .* exp(-1i*M/(1-M^2)*omega/c.*(x-xS));
souVx = M / (rho*c) * souU;

end