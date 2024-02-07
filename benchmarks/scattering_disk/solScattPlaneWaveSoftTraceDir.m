function [ Val ] = solScattPlaneWaveSoftTraceDir(k,R,thetaTab)

nEnd = floor(k*R) + 10;

% Compute solution

Val = zeros(size(thetaTab));
for n = 0:nEnd
    Val = Val - (1+(n>0)) * (1i)^n * besselj(n,k*R) * cos(n*thetaTab);
end

end