function [ Val ] = solScattModeHard(k,R,m,xTab,yTab)

zTab = xTab + 1i*yTab;
zTab = zTab.*(abs(zTab)>1);
rTab = abs(zTab);
thetaTab = angle(zTab);

% Compute Hankel functions

Hankel = zeros(1,m+2);
for n = 0:(m+1)
    Hankel(n+1) = besselj(n,k*R) + 1i * bessely(n,k*R);
end
dHankel = zeros(1,m+1);
dHankel(1) = -Hankel(2);
for n = 1:m
    dHankel(n+1) = Hankel(n) - n/(k*R) * Hankel(n+1);
end

Hankelr = besselj(m,k*rTab) + 1i * bessely(m,k*rTab);

% Compute solution

Val = real(dHankel(m+1))/dHankel(m+1) * Hankelr .* exp(1i*m*thetaTab);

end