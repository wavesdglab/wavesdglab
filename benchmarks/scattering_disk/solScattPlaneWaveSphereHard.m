function [ Val ] = solScattPlaneWaveSphereHard(k,R,xTab,yTab,zTab)

rTab = sqrt(xTab.^2+yTab.^2+zTab.^2);
nEnd = floor(k*R) + 20;

% Compute Legendre functions

kDir = [1 1 0]/sqrt(2);

varLeg = (kDir(1)*xTab+kDir(2)*yTab+kDir(3)*zTab)./rTab;
varLeg(rTab < 1) = NaN;
Legendre = zeros(size(varLeg,1),size(varLeg,2),size(varLeg,3),nEnd+1);
Legendre(:,:,:,1) = 1;
Legendre(:,:,:,2) = varLeg;
for n = 1:(nEnd-1)
    Legendre(:,:,:,n+2) = ((2*n+1) * varLeg .* Legendre(:,:,:,n+1) - n*Legendre(:,:,:,n))/(n+1);
end

% Compute derivative of spherical Hankel functions

Hankel = zeros(1,nEnd+2);
for n = 0:(nEnd+1)
    Hankel(n+1) = sqrt(pi/(2*k*R)) * (besselj(n+0.5,k*R) + 1i * bessely(n+0.5,k*R));
end
dHankel = zeros(1,nEnd+1);
dHankel(1) = -Hankel(2);
for n = 0:nEnd
    dHankel(n+1) = n/(k*R) * Hankel(n+1) - Hankel(n+2);
end

% Compute solution

Val = zeros(size(zTab));
for n = 0:nEnd
    
    % Spherical Hankel function (h_n^1(kr))
    HankelTab = sqrt(pi./(2*k*rTab)) .* (besselj(n+0.5,k*rTab) + 1i * bessely(n+0.5,k*rTab));
    
    coef = - (1i)^n * (2*n+1) * real(Hankel(n+1))/dHankel(n+1);
    Val = Val + coef * HankelTab .* Legendre(:,:,:,n+1);
end

end