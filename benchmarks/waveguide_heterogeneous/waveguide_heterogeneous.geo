// Gmsh project created on Fri Oct 13 14:07:11 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2, 0, 0, 1.0};
//+
Point(3) = {4, 0, 0, 1.0};
//+
Point(4) = {4, 1, 0, 1.0};
//+
Point(5) = {2, 1, 0, 1.0};
//+
Point(6) = {0, 1, 0, 1.0};
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
Physical Curve(1) = {1, 2};
//+
Physical Curve(2) = {3};
//+
Physical Curve(3) = {4, 5};
//+
Physical Curve(4) = {6};
//+
Physical Surface(12) = {1, 2};
