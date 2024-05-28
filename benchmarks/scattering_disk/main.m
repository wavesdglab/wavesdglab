close all;
clear all;

k = 20.;
R = 1;

Npt = 200;
d = 1.1;

X = linspace(-d,d,Npt);
Y = linspace(-d,d,Npt);
x = zeros(Npt,Npt);
y = zeros(Npt,Npt);
for i=1:Npt
    x(i,:) = X(i);
    y(:,i) = Y(i);
end
r = sqrt(x.^2+y.^2);

p = solScattPlaneWaveHard(k,R,x,y);
% p = solScattPlaneWaveSoft(k,R,x,y);

% mode = 20;
% p = solScattModeHard(k,R,mode,x,y);
% p = solScattModeSoft(k,R,mode,x,y);

p(r<R) = NaN;
pcolor(x,y,real(p));
shading interp;
colorbar;