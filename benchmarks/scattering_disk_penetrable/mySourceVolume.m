function val = mySourceVolume(x,y)

global LdomX LdomY Rdom omega cAir cObj Rdisk

val = zeros(size(x));
% if(~isempty(LdomX) && ~isempty(LdomY))  %%% For total formulation
%     if ((mean(abs(x)) >= LdomX) || (mean(abs(y)) >= LdomY))
%         kAir = omega/cAir;
%         val = exp(1i*kAir*x);
%     end
% end
% if(~isempty(Rdom))
%     r = sqrt(x.*x + y.*y);
%     if (mean(r) >= Rdom)
%         kAir = omega/cAir;
%         val = exp(1i*kAir*x);
%     end
% end

%%% For scattering formulation
r = sqrt(x.*x + y.*y);
if (mean(r) <= Rdisk)
    nhet = 1 /cObj;
    kAir = omega/cAir;
    val = exp(1i*kAir*x) .* omega^2 * (1 - nhet^2);
end

end