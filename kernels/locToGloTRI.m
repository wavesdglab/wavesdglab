function [x, y] = locToGloTRI(u, v, V1, V2, V3)

l1 =  (v+1)/2;
l2 = -(u+v)/2;
l3 =  (u+1)/2;
x = V1(1)*l1 + V2(1)*l2 + V3(1)*l3;
y = V1(2)*l1 + V2(2)*l2 + V3(2)*l3;

% dx/du = 0.5*(V3(1)-V2(1))
% dx/dv = 0.5*(V1(1)-V2(1))
% dy/du = 0.5*(V3(2)-V2(2))
% dy/dv = 0.5*(V1(2)-V2(2))
% J = [V3-V2 V1-v2] * 0.5

end