function [c, rho] = myCoef(x,y)

global cAir cObj rhoAir rhoObj Rdisk

c = cAir;
rho = rhoAir;
if(sqrt(x*x+y*y) < Rdisk)
    c = cObj;
    rho = rhoObj;
end

end