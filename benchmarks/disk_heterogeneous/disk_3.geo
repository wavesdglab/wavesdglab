//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, -0.5, 0, 1.0};
//+
Point(3) = {0, 0.5, 0, 1.0};
//+
Point(4) = {0, -0.25, 0, 1.0};
//+
Point(5) = {0, 0.25, 0, 1.0};
//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 0.25, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 0.5, 0, 2*Pi};
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
Field[1] = Distance;
Field[1].CurvesList = {2};
Field[1].Sampling = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.05;
Field[2].SizeMax = 0.06;
Field[2].DistMin = 0.2;
Field[2].DistMax = 0.2;
Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 5;