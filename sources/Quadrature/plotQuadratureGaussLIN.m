function plotQuadratureGaussLIN(degree)

[x,weights] = quadratureGaussLIN(degree);
scatter(x,zeros(size(x,1),size(x,2)),50,weights);
colorbar;

end