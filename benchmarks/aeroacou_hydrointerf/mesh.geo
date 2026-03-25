//Mesh.MshFileVersion = 2.2;

Point(1) = {-1, -1, 0};
Point(2) = {1, -1, 0};
Point(3) = {1, -0.000001, 0};
Point(4) = {0, 0, 0};
Point(5) = {1, 0.000001, 0};
Point(6) = {1, 1, 0};
Point(7) = {-1, 1, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Line(7) = {7};
Physical Surface(1) = {1};