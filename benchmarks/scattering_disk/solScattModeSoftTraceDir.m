function [ Val ] = solScattModeSoftTraceDir(k,R,m,thetaTab)

% Compute solution

Val = besselj(m,k*R) * exp(1i*m*thetaTab);

end