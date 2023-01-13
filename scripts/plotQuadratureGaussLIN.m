close all;
clear all;

Q = 10;
[nodes, weights] = quadratureGaussLIN(Q);

scatter(nodes,zeros(size(nodes,1),size(nodes,2)),50,weights);
colorbar;
