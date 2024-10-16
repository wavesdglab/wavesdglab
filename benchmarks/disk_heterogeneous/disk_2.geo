// Gmsh project created on Thu Feb 29 14:32:22 2024
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.5, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 1, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Surface(1) = {1};
//+
Curve Loop(3) = {2};
//+
Curve Loop(4) = {1};
//+
Surface(2) = {3, 4};
//+
Physical Curve(5) = {2};
//+
Physical Surface(6) = {1};
//+
Curve Loop(5) = {2};
//+
Curve Loop(6) = {1};
//+
Plane Surface(2) = {5, 6};
//+
Physical Surface(7) = {2};
