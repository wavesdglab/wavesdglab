close all;
clear all;

headers2D()

degree = 4;

x = -1:0.01:1;
y = -1:0.01:1;
[X,Y] = meshgrid(x,y);
for i=1:size(x,2)
    for j=1:size(y,2)
        if(x(i)+y(j)>0)
            X(i,j) = NaN;
            Y(i,j) = NaN;
        end
    end
end
val = functionsShapeTRI(X(:),Y(:),degree);
[valDX, valDY] = functionsShapeDerTRI(X(:),Y(:),degree);

N = 3 + 3*(degree-1) + (degree-1)*(degree-2)/2;
subplotX = ceil(N/3);
subplotY = 3;

figure(1);
title('Shape Functions');
for n=1:size(val,2)
    C = reshape(val(:,n),size(x,2),size(y,2));
    subplot(subplotX,subplotY,n);
    %h = pcolor(X,Y,C); set(h,'EdgeColor','none');
    %h = surf(X,Y,C); set(h,'EdgeColor','none');
    h = contourf(X,Y,C);
    colorbar;
    box off;
end

figure(2);
title('Shape Functions (derivative DX)');
for n=1:size(val,2)
    C = reshape(valDX(:,n),size(x,2),size(y,2));
    subplot(subplotX,subplotY,n);
    %h = pcolor(X,Y,C); set(h,'EdgeColor','none');
    %h = surf(X,Y,C); set(h,'EdgeColor','none');
    h = contourf(X,Y,C);
    colorbar;
    box off;
end

figure(3);
title('Shape Functions (derivative DY)');
for n=1:size(val,2)
    C = reshape(valDY(:,n),size(x,2),size(y,2));
    subplot(subplotX,subplotY,n);
    %h = pcolor(X,Y,C); set(h,'EdgeColor','none');
    %h = surf(X,Y,C); set(h,'EdgeColor','none');
    h = contourf(X,Y,C);
    colorbar;
    box off;
end