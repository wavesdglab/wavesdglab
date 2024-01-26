//Mesh.MshFileVersion = 2.2;

Printf("Parameters : LdomX = %f, LdomY = %f, LpmlX = %f, LpmlY = %f, Rdisk = %f\n", LdomX, LdomY, LpmlX, LpmlY, Rdisk);

// Points for exterior boundary
Point(1) = {-(LdomX+LpmlX), -(LdomY+LpmlY), 0};
Point(2) = { (LdomX+LpmlX), -(LdomY+LpmlY), 0};
Point(3) = { (LdomX+LpmlX),  (LdomY+LpmlY), 0};
Point(4) = {-(LdomX+LpmlX),  (LdomY+LpmlY), 0};

// Points for interior disk
Point(5) = {0, -Rdisk, 0};
Point(6) = { Rdisk, 0, 0};
Point(7) = {0,  Rdisk, 0};
Point(8) = {-Rdisk, 0, 0};

// Point for center of disk
Point(9) = {0, 0, 0};

// Points for PML/domain interface
Point(10) = {-LdomX, -LdomY, 0};
Point(11) = { LdomX, -LdomY, 0};
Point(12) = { LdomX,  LdomY, 0};
Point(13) = {-LdomX,  LdomY, 0};

// Exterior boundary
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Curve Loop(1) = {1, 2, 3, 4};

// Boundary of the disk
Circle(5) = {6, 9, 7};
Circle(6) = {7, 9, 8};
Circle(7) = {8, 9, 5};
Circle(8) = {5, 9, 6};
Curve Loop(2) = {5, 6, 7, 8};

// PML/domain interface
Line(9) = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 13};
Line(12) = {13, 10};
Curve Loop(3) = {11, 12, 9, 10};

// Definition of the surfaces
Plane Surface(1) = {3, 2};
Plane Surface(2) = {1, 3};

// Definition of the physical regions
Physical Curve("bd_pml", 1) = {1,2,3,4};
Physical Curve("bd_obstacle", 2) = {5,6,7,8};
Physical Surface("domain", 1) = {1};
Physical Surface("pml", 2) = {2};