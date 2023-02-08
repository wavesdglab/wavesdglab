% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function plotQuadratureGaussLIN(degree)

[x,weights] = quadratureGaussLIN(degree);
scatter(x,zeros(size(x,1),size(x,2)),50,weights);
colorbar;

end