Printf("Parameters: LdomX = %f, LdomY = %f, LpmlX = %f, LpmlY = %f, Rdisk = %f\n", LdomX, LdomY, LpmlX, LpmlY, Rdisk);

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

// Boundary of the disk
Circle(5) = {6, 9, 7};
Circle(6) = {7, 9, 8};
Circle(7) = {8, 9, 5};
Circle(8) = {5, 9, 6};

// PML/domain interface
Line(9) = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 13};
Line(12) = {13, 10};

// Domain (inside disk)
Curve Loop(1) = {5, 6, 7, 8};
Plane Surface(1) = {1};

// Domain (outside disk)
Curve Loop(2) = {5, 6, 7, 8, 11, 12, 9, 10};
Plane Surface(2) = {2};

// Layer
Curve Loop(3) = {-1, -2, -3, -4, 11, 12, 9, 10};
Plane Surface(3) = {3};

// Definition of the physical regions
Physical Curve(201) = {1,2,3,4};  // Exterior boundary of the layer
Physical Surface(301) = {1};      // Domain (inside disk)
Physical Surface(302) = {2};      // Domain (outside disk)
Physical Surface(303) = {3};      // Layer