function [ Val ] = modeTurbofan(k,R,m,xTab,yTab)

zTab = xTab + 1i*yTab;
zTab = zTab.*(abs(zTab)>1);
rTab = abs(zTab);
thetaTab = angle(zTab);

% Compute Hankel and dHankel functions for r = R
Hankel = zeros(1,m+2);
for n = 0:(m+1)
    Hankel(n+1) = besselj(n,k*R) + 1i * bessely(n,k*R);
end
dHankel = zeros(1,m+1);
dHankel(1) = -Hankel(2);
for n = 1:m
    dHankel(n+1) = Hankel(n) - n/(k*R) * Hankel(n+1);
end

% Compute Hankel functions for {r}

Hankelr = besselj(m,k*rTab) + 1i * bessely(m,k*rTab);

% Compute solution

Val = ( real(Hankelr) + real(dHankel(m+1))/imag(dHankel(m+1)) * imag(Hankelr) ) .* exp(1i*m*thetaTab);

end