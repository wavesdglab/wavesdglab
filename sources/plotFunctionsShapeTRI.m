% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function plotFunctionsShapeTRI()

degree = 4;

r = -1:0.01:1;
s = -1:0.01:1;
[R,S] = meshgrid(r,s);
for i=1:size(r,2)
    for j=1:size(s,2)
        if(r(i)+s(j)>0)
            R(i,j) = NaN;
            S(i,j) = NaN;
        end
    end
end
val = functionsShapeTRI(R(:),S(:),degree);
[valDX, valDY] = functionsShapeDerTRI(R(:),S(:),degree);

N = 3 + 3*(degree-1) + (degree-1)*(degree-2)/2;
subplotX = ceil(N/3);
subplotY = 3;

figure(3);
title('Shape Functions');
for n=1:size(val,2)
    C = reshape(val(:,n),size(r,2),size(s,2));
    subplot(subplotX,subplotY,n);
    %h = pcolor(R,S,C); set(h,'EdgeColor','none');
    %h = surf(R,S,C); set(h,'EdgeColor','none');
    h = contourf(R,S,C);
    colorbar;
    box off;
end

% figure(2);
% title('Shape Functions (derivative DX)');
% for n=1:size(val,2)
%     C = reshape(valDX(:,n),size(r,2),size(s,2));
%     subplot(subplotX,subplotY,n);
%     %h = pcolor(R,S,C); set(h,'EdgeColor','none');
%     %h = surf(R,S,C); set(h,'EdgeColor','none');
%     h = contourf(R,S,C);
%     colorbar;
%     box off;
% end
% 
% figure(3);
% title('Shape Functions (derivative DY)');
% for n=1:size(val,2)
%     C = reshape(valDY(:,n),size(r,2),size(s,2));
%     subplot(subplotX,subplotY,n);
%     %h = pcolor(R,S,C); set(h,'EdgeColor','none');
%     %h = surf(R,S,C); set(h,'EdgeColor','none');
%     h = contourf(R,S,C);
%     colorbar;
%     box off;
% end

end