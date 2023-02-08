% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [uQ, vQ, weights] = quadratureGaussTRI(degree)

[nodes, weights] = triasymq(degree, [-1 -1], [1 -1], [-1 1]);

uQ = nodes(1,:)';
vQ = nodes(2,:)';
weights = weights';

end