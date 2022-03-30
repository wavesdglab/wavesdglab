global L;
global k;
global theta;

L     = 4.;
k     = 4.*pi;
theta = 10.*(pi/180.);

N = 128;

x = linspace(0,L,L*N);
y = linspace(0,1,  N);

[X,Y] = meshgrid(x,y);

[u,dx,dy] = waveguide(X,Y);

pcolor(X,Y,real(u));
