// Gmsh project created on Thu Sep 28 11:22:50 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, -0, 1.0};
//+
Point(2) = {0.5, 0, -0, 1.0};
//+
Point(3) = {1, 0, -0, 1.0};
//+
Point(4) = {1, 1, -0, 1.0};
//+
Point(5) = {0.5, 1, -0, 1.0};
//+
Point(6) = {0, 1, -0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Line(7) = {2, 5};
//+
Curve Loop(1) = {1, 7, 5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 3, 4, -7};
//+
Plane Surface(2) = {2};
//+
Physical Curve(4) = {1, 2};
//+
Physical Curve(1) = {3};
//+
Physical Curve(2) = {4, 5};
//+
Physical Curve(3) = {6};
//+
Physical Surface(12) = {1, 2};
