function [ Val ] = solScattPlaneWaveSoft(k,R,xTab,yTab)

zTab = xTab + 1i*yTab;
zTab = zTab.*(abs(zTab)>0.5);
rTab = abs(zTab);
thetaTab = angle(zTab);

nEnd = floor(k*R) + 10;

% Compute Hankel functions

Hankel = zeros(1,nEnd+2);
for n = 0:nEnd
    Hankel(n+1) = besselj(n,k*R) + 1i * bessely(n,k*R);
end

% Compute solution

Val = zeros(size(zTab));
for n = 0:nEnd
    tmp = - (1+(n>0)) * (1i)^n * real(Hankel(n+1))/Hankel(n+1);
    HankelTab = besselj(n,k*rTab) + 1i * bessely(n,k*rTab);
    Val = Val + tmp * HankelTab .* cos(n*thetaTab);
end

end