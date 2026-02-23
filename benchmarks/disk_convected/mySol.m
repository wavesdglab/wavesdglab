function [solU, solDx, solDy, solVx, solVy] = mySol(x,y)

global v0d rho c omega;

xs = -v0d;
ys = 0;
xg = x-xs;
yg = y-ys;

alpha = omega/(c^2-v0d^2);

coef = 1i / (4*c*sqrt(c^2-v0d^2)) .* exp(-1i*alpha*v0d .* xg);
dist = sqrt(c^2*xg.^2 + (c^2-v0d^2).*yg.^2);
D = alpha*dist;
dDdx = alpha*c^2*xg./dist;
dDdy = alpha*(c^2-v0d^2)*yg./dist;
G = coef .* besselh(0,D);
dGdx = -coef .* (dDdx .* besselh(1,D) + 1i*alpha*v0d .* besselh(0,D));
dGdy = -coef .* (dDdy .* besselh(1,D));

solU = rho*c^2 * (-1i*omega*G + v0d*c*dGdx);
solVx = - c^2 * dGdx;
solVy = - c^2 * dGdy;
solDx = zeros(size(x));
solDy = zeros(size(x));

end