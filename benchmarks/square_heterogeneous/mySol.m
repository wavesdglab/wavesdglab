function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global k1 k2 k3 eta1 eta2 eta3

% Coordinates of the center of the narrow Gaussian point source
x0=0;
y0=0;

X=mean(x);
Y=mean(y);

if Y>0.25
    k=k1;
    eta=eta1;
end
if Y>0 && Y<0.25
    k=k2;
    eta=eta2;
end
if Y>-0.25 && Y<0
    if X>-0.25 && X<0.25
        k=k2;
        eta=eta2;
    else
        k=k3;
        eta=eta3;
    end
end
if Y<-0.25
    k=k3;
    eta=eta3;
end
if Y==0
    k=k2;
    eta=eta2;
end
if Y==-0.25
    k=k3;
    eta=eta3;
end

solU = 0*x + 0*y;
solDx = 0*x + 0*y;
solDy = 0*x + 0*y;
solF = eta*exp(-(4*k/pi)^2*((x-x0).^2+(y-y0).^2));
solVx = 1/(1i*k*eta) * solDx;
solVy = 1/(1i*k*eta) * solDy;

end