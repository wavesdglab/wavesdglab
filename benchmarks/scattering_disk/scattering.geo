//Mesh.MshFileVersion = 2.2;

Point(1) = {-1, -1, 0};
Point(2) = {1, -1, 0};
Point(3) = {1, 1, 0};
Point(4) = {-1, 1, 0};
Point(5) = {0, -0.5, 0};
Point(6) = {0.5, 0, 0};
Point(7) = {0, 0.5, 0};
Point(8) = {-0.5, 0, 0};
Point(9) = {0, 0, 0};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line Loop(1) = {1, 2, 3, 4};

Circle(5) = {6, 9, 7};
Circle(6) = {7, 9, 8};
Circle(7) = {8, 9, 5};
Circle(8) = {5, 9, 6};

Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1, 2};
Physical Curve(1) = {1,2,3,4};
Physical Curve(2) = {5,6,7,8};
Physical Surface(1) = {1};