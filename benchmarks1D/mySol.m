function val = mySol(x)
global k BCLeft BCRight

if(strcmp(BCLeft,'PER') && strcmp(BCRight,'PER'))
    val = exp(1i*20*pi*x);
elseif(strcmp(BCLeft,'DIR') && strcmp(BCRight,'DIR'))
    val = 2*sin(k/2)/sin(k) * cos(k*(x-1/2)) - 1;
elseif(strcmp(BCLeft,'DIR') && strcmp(BCRight,'ABC'))
    val = exp(1i*k*x) - 1i*exp(1i*k)*sin(k*x) - 1;
elseif(strcmp(BCLeft,'ABC') && strcmp(BCRight,'ABC'))
    val = exp(1i*k/2)*cos(k*(x-1/2)) - 1;
elseif(strcmp(BCLeft,'NEU') && strcmp(BCRight,'ABC'))
    val = exp(1i*k*x);
else
    val = exp(1i*k/2)*cos(k*(x-1/2)) - 1;
end

end