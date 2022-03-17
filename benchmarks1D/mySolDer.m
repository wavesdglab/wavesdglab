function val = mySolDer(x)
global k BCLeft BCRight

if(strcmp(BCLeft,'PER') && strcmp(BCRight,'PER'))
    val = 1i*20*pi* exp(1i*20*pi*x);
elseif(strcmp(BCLeft,'DIR') && strcmp(BCRight,'DIR'))
    val = -k * 2*sin(k/2)/sin(k) * sin(k*(x-1/2));
elseif(strcmp(BCLeft,'DIR') && strcmp(BCRight,'ABC'))
    val = 1i*k * (exp(1i*k*x) - exp(1i*k)*cos(k*x));
elseif(strcmp(BCLeft,'ABC') && strcmp(BCRight,'ABC'))
    val = -k * exp(1i*k/2) * sin(k*(x-1/2));
elseif(strcmp(BCLeft,'NEU') && strcmp(BCRight,'ABC'))
    val = 1i*k*exp(1i*k*x);
else
    val = -k * exp(1i*k/2)*sin(k*(x-1/2));
end

end