function [solP, derP, solV, sou, souP, souV] = mySol1D_heterogeneous(x,k,eta)
global omega BCLeft BCRight

if(strcmp(BCLeft,'PER') && strcmp(BCRight,'PER'))
    solP = exp(1i*20*pi*x);
    derP = 1i*20*pi* exp(1i*20*pi*x);
    sou  = ((20*pi)^2-k^2) * exp(1i*20*pi*x);
elseif(strcmp(BCLeft,'DIR') && strcmp(BCRight,'DIR'))
%     solP = 2*sin(k/2)/sin(k) * cos(k*(x-1/2)) - 1;
%     derP = -k * 2*sin(k/2)/sin(k) * sin(k*(x-1/2));
%     sou  = k^2*ones(size(x,1),size(x,2));

%     c = omega / k;
%     solP = exp(1i*2*pi*x)-1;
%     derP = 1i*2*pi*exp(1i*2*pi*x);
%     sou =  k^2 * (solP * (c^2 - 1) + c^2);

%     c = omega / k;
%     solP = eta * (exp(1i*4*pi*x)-1);
%     derP = eta * 1i * 4 * pi * exp(1i*4*pi*x);
%     sou = 16 * pi^2 * eta * (exp(1i*4*pi*x)-1/c^2*exp(1i*4*pi*x)+1/c^2);

%     c = omega / k;
%     rho = eta / c;
%     solP = rho * (exp(1i*4*pi*x)-1);
%     derP = rho * 1i * 4 * pi * exp(1i*4*pi*x);
%     sou = 16 * pi^2 * rho * (exp(1i*4*pi*x)-1/c^2*exp(1i*4*pi*x)+1/c^2);

% The solution is continuous, but it is defined piecewise
rho1 = 1;
c1 = 2;

rho2 = 1;
c2 = 3;

k1 = omega / c1;
k2 = omega / c2;

eta1 = rho1 * c1;
eta2 = rho2 * c2;

I = 0.5;              % amplitude of the incident pressure field
if (min(x)<1/2)
    solP = I*(exp(1i*k1*x)+(eta2-eta1)/(eta1+eta2)*exp(1i*k1)*exp(-1i*k1*x));
    derP = 1i*k1*I*(exp(1i*k1*x)-(eta2-eta1)/(eta1+eta2)*exp(1i*k1)*exp(-1i*k1*x));
    sou = 0*x;
else
    solP = I*(2*eta2/(eta1+eta2)*exp(1i*(k1-k2)/2)*exp(1i*k2*x));
    derP = 1i*k2*I*2*eta2/(eta1+eta2)*exp(1i*(k1-k2)/2)*exp(1i*k2*x);
    sou = 0*x;
end


elseif(strcmp(BCLeft,'DIR') && strcmp(BCRight,'ABC'))
%     solP = exp(1i*k*x) - 1i*exp(1i*k)*sin(k*x) - 1;
%     derP = 1i*k * (exp(1i*k*x) - exp(1i*k)*cos(k*x));
%     sou  = k^2*ones(size(x,1),size(x,2));

% The solution is continuous, but it is defined piecewise
rho1 = 1;
c1 = 0.1;

rho2 = 1;
c2 = 2;

k1 = omega / c1;
k2 = omega / c2;

eta1 = rho1 * c1;
eta2 = rho2 * c2;

I = 0.5;              % amplitude of the incident pressure field
if (min(x)<1/2)
    solP = I*(exp(1i*k1*x)+(eta2-eta1)/(eta1+eta2)*exp(1i*k1)*exp(-1i*k1*x));
    derP = 1i*k1*I*(exp(1i*k1*x)-(eta2-eta1)/(eta1+eta2)*exp(1i*k1)*exp(-1i*k1*x));
    sou = 0*x;
else
    solP = I*(2*eta2/(eta1+eta2)*exp(1i*(k1-k2)/2)*exp(1i*k2*x));
    derP = 1i*k2*I*2*eta2/(eta1+eta2)*exp(1i*(k1-k2)/2)*exp(1i*k2*x);
    sou = 0*x;
end

elseif(strcmp(BCLeft,'ABC') && strcmp(BCRight,'ABC'))
%     solP = exp(1i*k/2)*cos(k*(x-1/2)) - 1;
%     derP = -k * exp(1i*k/2) * sin(k*(x-1/2));
%     sou  = k^2*ones(size(x,1),size(x,2));

%     solP = exp(1i*k*x)-1/(1i*k);
%     derP = 1i*k*exp(1i*k*x);
%     sou = 0*x;

% The solution is continuous, but it is defined piecewise
rho1 = 1;
c1 = 0.1;

rho2 = 1;
c2 = 2;

k1 = omega / c1;
k2 = omega / c2;

eta1 = rho1 * c1;
eta2 = rho2 * c2;

I = 0.5;              % amplitude of the incident pressure field
if (min(x)<1/2)
    solP = I*(exp(1i*k1*x)+(eta2-eta1)/(eta1+eta2)*exp(1i*k1)*exp(-1i*k1*x));
    derP = 1i*k1*I*(exp(1i*k1*x)-(eta2-eta1)/(eta1+eta2)*exp(1i*k1)*exp(-1i*k1*x));
    sou = 0*x;
else
    solP = I*(2*eta2/(eta1+eta2)*exp(1i*(k1-k2)/2)*exp(1i*k2*x));
    derP = 1i*k2*I*2*eta2/(eta1+eta2)*exp(1i*(k1-k2)/2)*exp(1i*k2*x);
    sou = 0*x;
end

elseif(strcmp(BCLeft,'NEU') && strcmp(BCRight,'ABC'))
    solP = exp(1i*k*x);
    derP = 1i*k*exp(1i*k*x);
    sou  = k^2*ones(size(x,1),size(x,2));
else
    solP = exp(1i*k/2)*cos(k*(x-1/2)) - 1;
    derP = -k * exp(1i*k/2)*sin(k*(x-1/2));
    sou  = k^2*ones(size(x,1),size(x,2));
end

solV = 1/(1i*k*eta)*derP;
souP = -1/(1i*k*eta)*sou;
souV = zeros(size(x,1),size(x,2));

end