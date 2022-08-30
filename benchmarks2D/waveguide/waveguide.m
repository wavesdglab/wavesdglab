function [u,dx,dy] = waveguide(x,y,k,L,theta)
N = 48;

u  = zeros(size(x));
dx = zeros(size(x));
dy = zeros(size(x));

for n=1:N
    npi = n*pi;
    g = waveguide_planewave_coefficient(n,k,theta);
    [v,d] = waveguide_ode(n,g,x,k,L);
    u  = u  +     v.*sin(npi*y);
    dx = dx +     d.*sin(npi*y);
    dy = dy + npi*v.*cos(npi*y);
end
end