x = linspace(0,1,50);
y = x';
m = 2; n = 2;
% k = 2.01*sqrt(2)*pi;
k = 1;
z = (sin(n*x*pi).*sin(m*y*pi) + sin(m*x*pi).*sin(n*y*pi))./2;%.*  16*m*n*pi^2 /(pi^2*(m^2 + n^2)-k^2)
[X,Y] = meshgrid(y,x);
surf(X,Y,z)