function [ Val ] = scattDiskPenetrable_Solution(kAir,kObj,rhoAir,rhoObj,R,xTab,yTab)

zTab = xTab + 1i*yTab;
rTab = abs(zTab);
tTab = angle(zTab);

%nEnd = floor(max(kAir,kObj)*R) + 10;
nEnd = 30;

% Compute functions

HankelAir = zeros(1,nEnd+2);
BesselObj = zeros(1,nEnd+2);
for n = 0:(nEnd+1)
    HankelAir(n+1) = besselj(n,kAir*R) + 1i * bessely(n,kAir*R);
    BesselObj(n+1) = besselj(n,kObj*R);
end
dHankelAir = zeros(1,nEnd+1);
dBesselObj = zeros(1,nEnd+1);
dHankelAir(1) = -HankelAir(2);
dBesselObj(1) = -BesselObj(2);
for n = 1:nEnd
    dHankelAir(n+1) = HankelAir(n) - n/(kAir*R) * HankelAir(n+1);
    dBesselObj(n+1) = BesselObj(n) - n/(kObj*R) * BesselObj(n+1);
end

% Compute solution

ValAir = zeros(size(zTab));
ValObj = zeros(size(zTab));
alpha = (rhoAir*kObj)/(rhoObj*kAir);
for n = 0:nEnd
    tmp = - (1+(n>0)) * (1i)^n;
    numAir = real(dHankelAir(n+1)) * BesselObj(n+1) - alpha * real(HankelAir(n+1)) * dBesselObj(n+1);
    numObj = real(dHankelAir(n+1)) * HankelAir(n+1) -         real(HankelAir(n+1)) * dHankelAir(n+1);
    denum = dHankelAir(n+1) * BesselObj(n+1) - alpha * HankelAir(n+1) * dBesselObj(n+1);
    HankelAirTab = besselj(n,kAir*rTab) + 1i * bessely(n,kAir*rTab);
    BesselObjTab = besselj(n,kObj*rTab);
    ValAir = ValAir + tmp * (numAir/denum) * HankelAirTab .* cos(n*tTab);
    ValObj = ValObj + tmp * (numObj/denum) * BesselObjTab .* cos(n*tTab);
end

% Assemble solution

Val = zeros(size(zTab));
for i=1:size(xTab,1)
    for j=1:size(xTab,2)
        ValInc = exp(1i*kAir*xTab(i,j));
        if(rTab(i,j) < R)
            Val(i,j) = ValObj(i,j);
        else
            Val(i,j) = ValAir(i,j) + ValInc;
        end
    end
end