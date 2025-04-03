function [souU, souDx, souDy, souVx, souVy] = mySourceSurface(x,y)

global M c omega;

souU = zeros(size(x));
souDx = zeros(size(x));
souDy = zeros(size(x));
souVx = zeros(size(x));
souVy = zeros(size(x));

A = 1; % amplitude of the vorticity wave
H = 1; % width of the waveguide
n = 1; % mode

souU = 0.*x +0.*y;
souVx = souVx + A * (1i*n*pi*M*c)/(omega*H) * cos(n*pi*y/H) .* exp(1i*omega*x/(M*c));
souVy = souVy + A * sin(n*pi*y/H) .* exp(1i*omega*x/(M*c));