SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
Circle(2) = {-0, -0, 0, 1, 0, 2*Pi};
Curve Loop(1) = {1};
Plane Surface(1) = {1};
Physical Curve(3) = {1};
Physical Surface(4) = {1};