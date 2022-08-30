function [ Val ] = solScattModeHardTraceNeu(k,R,m,thetaTab)

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

% Compute solution

Val = k * real(dHankel(m+1)) * exp(1i*m*thetaTab);

end