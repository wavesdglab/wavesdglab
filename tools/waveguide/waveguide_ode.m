function [v,d] = waveguide_ode(n,g,x,k,L)
threshold = 1.e-10;

kappa = n*n*pi*pi-k*k;
ik    = 1i*k;

if (abs(kappa) > threshold)
    s = sqrt(kappa);
    
    rvL =   exp(s*L)-  exp(-s*L);
    rdL = s*exp(s*L)+s*exp(-s*L);
    
    rvx =   exp(s*x)-  exp(-s*x);
    rdx = s*exp(s*x)+s*exp(-s*x);
else
    rvL = L ;
    rdL = 1.;
    
    rvx = x ;
    rdx = 1.;
end

c = g/(rdL-ik*rvL);

v = c*rvx;
d = c*rdx;
end