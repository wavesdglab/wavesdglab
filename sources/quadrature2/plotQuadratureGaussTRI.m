function plotQuadratureGaussTRI(degree)

[x,y,weights] = quadratureGaussTRI(degree);
scatter(x,y,50,weights);
colorbar;

end