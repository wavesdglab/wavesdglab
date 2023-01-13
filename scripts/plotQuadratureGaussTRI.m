%close all;
clear all;

degree = 20;
[x, y, weights] = quadratureGaussTRI(degree);

scatter(x,y,50,weights);
colorbar;
