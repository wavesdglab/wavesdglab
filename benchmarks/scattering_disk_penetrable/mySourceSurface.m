function [souU, souDx, souDy, souVx, souVy] = mySourceSurface(x,~)

global omega cAir rhoAir
kAir = omega/cAir;

souU = exp(1i*kAir*x);
souDx = 1i*kAir*exp(1i*kAir*x);
souDy = 0;
souVx = souDx/(1i*omega*rhoAir);
souVy = 0;

end