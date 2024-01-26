function [ Val ] = solScattModeSoftTraceNeu(k,R,m,thetaTab)

% Compute Hankel functions

Hankel = zeros(1,m+2);
for n = 0:(m+1)
    Hankel(n+1) = besselj(n,k*R) + 1i * bessely(n,k*R);
end
if (m==0)
    dHankel = -Hankel(2);
else
    dHankel = Hankel(m) - m/(k*R) * Hankel(m+1);
end

% Compute solution

Val = k * real(Hankel(m+1))/Hankel(m+1) * dHankel * exp(1i*m*thetaTab);

end