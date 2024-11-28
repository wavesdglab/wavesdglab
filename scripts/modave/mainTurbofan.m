%close all;
clear all;

m = 4;
%m = 24;

Rin = 0.3586;
Rout = 1.2;

kr = 0:0.01:50;
dJout = besseljDer(m,kr*Rout);
dYout = besselyDer(m,kr*Rout);
dJin = besseljDer(m,kr*Rin);
dYin = besselyDer(m,kr*Rin);
VAL = dJout.*dYin - dJin.*dYout;
figure(1);
plot(kr,VAL);
ylim([-10,10])
krVal = [];
n = 0;
for i=1:size(VAL,2)-1
    if(VAL(i)*VAL(i+1) < 0)
        n = n+1;
        krVal(n) = kr(i);
    end
end

krVal'

kr = krVal(1);

Nvizu = 200;
x = -Rout:2/Nvizu:Rout;
y = -Rout:2/Nvizu:Rout;
[xTab,yTab] = meshgrid(x,y);
zTab = xTab + 1i*yTab;
zTab(abs(zTab)>Rout) = NaN;
zTab(abs(zTab)<Rin) = NaN;
rTab = abs(zTab);
thetaTab = angle(zTab);
%imagesc(x,y,thetaTab);

val = besselj(m,kr*rTab) - besseljDer(m,kr*Rin)/besselyDer(m,kr*Rin) * bessely(m,kr*rTab);
max(max(abs(val)))
val = val .* cos(m*thetaTab);
val = val/max(max(abs(val)));
figure(2);
imagesc(x,y,val);
colorbar;

Ncolor = 100;
map = [1 0 0
    1 1 1
    0 0 1.0];
map = [ones(Ncolor,1) (0:1/(Ncolor-1):1)' (0:1/(Ncolor-1):1)' ; (1:-1/(Ncolor-1):0)' (1:-1/(Ncolor-1):0)' ones(Ncolor,1)];
colormap(map);

% ---------------------

function val = besseljDer(nu,z)
if(nu == 0)
    val = -besselj(1,r);
else
    val = besselj(nu-1,z) - nu./z .* besselj(nu,z);
end
end

function val = besselyDer(nu,z)
if(nu == 0)
    val = -bessely(1,r);
else
    val = bessely(nu-1,z) - nu./z .* bessely(nu,z);
end
end