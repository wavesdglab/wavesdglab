%close all;
clear all;

headers2D()

degree = 20;
[x, y, weights] = quadratureGaussTRI(degree);

scatter(x,y,50,weights);
colorbar;
