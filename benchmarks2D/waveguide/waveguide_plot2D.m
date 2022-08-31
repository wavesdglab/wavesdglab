L     = 1.;
k     = 20.*pi;
theta = 18.*(pi/180.);

N = 128;

x = linspace(0,L,L*N);
y = linspace(0,1,  N);

[X,Y] = meshgrid(x,y);

[u,dx,dy] = waveguide(X,Y,k,L,theta);

pcolor(X,Y,real(u));
