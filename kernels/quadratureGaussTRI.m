function [x, y, weights] = quadratureGaussTRI(degree)

[rnodes, weights] = triasymq(degree, [-1 -1], [1 -1], [-1 1]);

x = rnodes(1,:)';
y = rnodes(2,:)';
weights = weights';

end