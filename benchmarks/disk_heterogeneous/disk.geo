//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, -1, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {0, -0.5, 0, 1.0};
//+
Point(5) = {0, 0.5, 0, 1.0};
//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 0.5, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 1, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2};
//+
Curve Loop(3) = {1};
//+
Plane Surface(2) = {2, 3};
//+
Physical Curve(4) = {2};
//+
Physical Surface(5) = {1};
//+
Physical Surface(6) = {2};
