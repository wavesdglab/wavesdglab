//Mesh.MshFileVersion = 2.2;

Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};
Line(1) = {1, 4}; // West
Line(2) = {4, 3}; // North
Line(3) = {3, 2}; // East
Line(4) = {2, 1}; // South
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Surface(1) = {1};