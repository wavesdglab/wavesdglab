% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function plotQuadratureGaussTRI(degree)

[x,y,weights] = quadratureGaussTRI(degree);
scatter(x,y,50,weights);
colorbar;

end