function [souU, souDx, souDy, souVx, souVy] = mySourceSurface(x,y)

global M rho c omega;

xs=-M;
ys=0;
k=omega/c;

souU  = zeros(size(x));
souDx = zeros(size(x));
souDy = zeros(size(x));
souVx = zeros(size(x));
souVy = zeros(size(x));

% A = 1i / (4*c^2*sqrt(1-M^2));
% z = k*sqrt(x.^2+(1-M^2).*y.^2)/(1-M^2);
% B = exp(-1i*M/(1-M^2)*k.*x);
% G = A .* besselh(0,z) .* B;
% dGdx = -A*k/(1-M^2) .* B .* ( besselh(1,z).*x./(sqrt(x.^2+(1-M^2).*y.^2)) + 1i*M*besselh(0,z));
% dGdy = -A*k.*y./(sqrt(x.^2+(1-M^2).*y.^2)) .* besselh(1,z) .* B;

A = 1i / (4*c^2*sqrt(1-M^2));
z = k*sqrt((x-xs).^2+(1-M^2).*(y-ys).^2)/(1-M^2);
B = exp(-1i*M/(1-M^2)*k.*(x-xs));
G = A .* besselh(0,z) .* B;
dGdx = -A*k/(1-M^2) .* B .* ( besselh(1,z).*(x-xs)./(sqrt((x-xs).^2+(1-M^2).*(y-ys).^2)) + 1i*M*besselh(0,z));
dGdy = -A*k.*(y-ys)./(sqrt((x-xs).^2+(1-M^2).*(y-ys).^2)) .* besselh(1,z) .* B;

souU = -1i*omega*G + M*c*dGdx;
souVx = - dGdx / rho;
souVy = - dGdy / rho;

end