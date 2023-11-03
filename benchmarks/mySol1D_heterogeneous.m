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

    c = omega / k;
    solP = eta * (exp(1i*4*pi*x)-1);
    derP = eta * 1i * 4 * pi * exp(1i*4*pi*x);
    sou = 16 * pi^2 * eta * (exp(1i*4*pi*x)-1/c^2*exp(1i*4*pi*x)+1/c^2);

elseif(strcmp(BCLeft,'DIR') && strcmp(BCRight,'ABC'))
    solP = exp(1i*k*x) - 1i*exp(1i*k)*sin(k*x) - 1;
    derP = 1i*k * (exp(1i*k*x) - exp(1i*k)*cos(k*x));
    sou  = k^2*ones(size(x,1),size(x,2));
elseif(strcmp(BCLeft,'ABC') && strcmp(BCRight,'ABC'))
    solP = exp(1i*k/2)*cos(k*(x-1/2)) - 1;
    derP = -k * exp(1i*k/2) * sin(k*(x-1/2));
    sou  = k^2*ones(size(x,1),size(x,2));
%     solP = exp(1i*k*x)-1/(1i*k);
%     derP = 1i*k*exp(1i*k*x);
%     sou = 0*x;
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