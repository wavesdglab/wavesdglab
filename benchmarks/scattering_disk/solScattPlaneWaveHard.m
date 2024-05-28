function [ Val ] = solScattPlaneWaveHard(k,R,xTab,yTab)

zTab = xTab + 1i*yTab;
rTab = max(abs(zTab),R);
thetaTab = angle(zTab);

nEnd = floor(k*R) + 10;

% Compute Hankel functions

Hankel = zeros(1,nEnd+2);
for n = 0:(nEnd+1)
    Hankel(n+1) = besselj(n,k*R) + 1i * bessely(n,k*R);
end
dHankel = zeros(1,nEnd+1);
dHankel(1) = -Hankel(2);
for n = 1:nEnd
    dHankel(n+1) = Hankel(n) - n/(k*R) * Hankel(n+1);
end

% Compute solution

Val = zeros(size(zTab));
for n = 0:nEnd
    tmp = - (1+(n>0)) * (1i)^n * real(dHankel(n+1))/dHankel(n+1);
    HankelTab = besselj(n,k*rTab) + 1i * bessely(n,k*rTab);
    Val = Val + tmp * HankelTab .* cos(n*thetaTab);
end

end