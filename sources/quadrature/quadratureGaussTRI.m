function [uQ, vQ, weights] = quadratureGaussTRI(degree)

[nodes, weights] = triasymq(degree, [-1 -1], [1 -1], [-1 1]);

uQ = nodes(1,:)';
vQ = nodes(2,:)';
weights = weights';

end