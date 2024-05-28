function [ Val ] = solScattModeSoft(k,R,m,xTab,yTab)

zTab = xTab + 1i*yTab;
zTab = zTab.*(abs(zTab)>1);
rTab = abs(zTab);
thetaTab = angle(zTab);

% Compute Hankel functions

HankelR = besselj(m,k*R) + 1i * bessely(m,k*R);
Hankelr = besselj(m,k*rTab) + 1i * bessely(m,k*rTab);

% Compute solution

Val = real(HankelR)/HankelR * Hankelr .* exp(1i*m*thetaTab);

end