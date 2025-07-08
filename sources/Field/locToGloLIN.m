function [x, y] = locToGloLIN(u, V1, V2)

l1 = (1-u)/2;
l2 = (1+u)/2;
x = V1(1)*l1 + V2(1)*l2;
y = V1(2)*l1 + V2(2)*l2;

% dx/du = 0.5*(V2(1)-V1(1))
% dy/du = 0.5*(V2(2)-V1(2))
% J = norm(V2-V1) * 0.5


end
