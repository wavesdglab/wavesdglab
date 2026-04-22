%close all;
clear all;

kAir = 5*pi;
kObj = 2.5*pi;
rhoAir = 1.;
rhoObj = 1.;
R = 1;

Npt = 200;
d = 1.5;
X = linspace(-d,d,Npt);
Y = linspace(-d,d,Npt);
x = zeros(Npt,Npt);
y = zeros(Npt,Npt);
for i=1:Npt
    x(i,:) = X(i);
    y(:,i) = Y(i);
end
r = sqrt(x.^2+y.^2);

p = scattDiskPenetrable_Solution(kAir,kObj,rhoAir,rhoObj,R,x,y,'scattered');

figure(1);
pcolor(x,y,real(p));
shading interp;
colorbar;