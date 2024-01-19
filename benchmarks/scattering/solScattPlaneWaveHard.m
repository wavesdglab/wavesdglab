function [ Val ] = solScattPlaneWaveHard(k,R_disk,xTab,yTab)

global L L_PML

zTab = xTab + 1i*yTab;
zTab = zTab.*(abs(zTab)>R_disk);
rTab = abs(zTab);
thetaTab = angle(zTab);

nEnd = floor(k*R_disk) + 10;

% Compute Hankel functions

Hankel = zeros(1,nEnd+2);
for n = 0:(nEnd+1)
    Hankel(n+1) = besselj(n,k*R_disk) + 1i * bessely(n,k*R_disk);
end
dHankel = zeros(1,nEnd+1);
dHankel(1) = -Hankel(2);
for n = 1:nEnd
    dHankel(n+1) = Hankel(n) - n/(k*R_disk) * Hankel(n+1);
end

% Compute solution

Val = zeros(size(zTab));
for n = 0:nEnd
    tmp = - (1+(n>0)) * (1i)^n * real(dHankel(n+1))/dHankel(n+1);
    HankelTab = besselj(n,k*rTab) + 1i * bessely(n,k*rTab);
    Val = Val + tmp * HankelTab .* cos(n*thetaTab);
end

end