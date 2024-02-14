function val = mySourceVolume(x,y)

global LdomX LdomY Rdom omega cAir

val = zeros(size(x));
if(~isempty(LdomX) && ~isempty(LdomY))
    if ((mean(abs(x)) >= LdomX) || (mean(abs(y)) >= LdomY))
        kAir = omega/cAir;
        val = exp(1i*kAir*x);
    end
end
if(~isempty(Rdom))
    r = sqrt(x.*x + y.*y);
    if (mean(r) >= Rdom)
        kAir = omega/cAir;
        val = exp(1i*kAir*x);
    end
end

end