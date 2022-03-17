function val = mySou(x)
global k BCLeft BCRight

if(strcmp(BCLeft,'PER') && strcmp(BCRight,'PER'))
    val = ((20*pi)^2-k^2) * exp(1i*20*pi*x);
else
    val = k^2*ones(size(x,1),size(x,2));
end

end