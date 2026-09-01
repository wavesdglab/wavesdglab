function [souU, souDx, souDy, souVx, souVy] = mySourceSurface(x,y)

global omega c rho v0;

souU = zeros(size(x));
souDx = zeros(size(x));
souDy = zeros(size(x));


lambda = 2*pi*v0(1)/omega;
sigma = lambda/sqrt(2);

gFun = exp(-y.^2/sigma^2);
fFun = 2*v0(1)*y/(1i*omega*sigma^2) .* gFun;

coef = rho*c*sqrt(v0(1));

souVx = coef * fFun .* exp(1i*omega*x/v0(1));
souVy = coef * gFun .* exp(1i*omega*x/v0(1));