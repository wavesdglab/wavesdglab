Printf("Parameters: Rdisk = %f, Rdom = %f, Rpml = %f\n", Rdisk, Rdom, Rpml);

Point(0) = {0, 0, 0};

// Boundary of scattering disk
Point(1) = {0, -Rdisk, 0};
Point(2) = { Rdisk, 0, 0};
Point(3) = {0,  Rdisk, 0};
Point(4) = {-Rdisk, 0, 0};
Circle(1) = {1, 0, 4};
Circle(2) = {4, 0, 3};
Circle(3) = {3, 0, 2};
Circle(4) = {2, 0, 1};

// Exterior boundary of the domain
Point(11) = {0, -Rdom, 0};
Point(12) = { Rdom, 0, 0};
Point(13) = {0,  Rdom, 0};
Point(14) = {-Rdom, 0, 0};
Circle(11) = {11, 0, 14};
Circle(12) = {14, 0, 13};
Circle(13) = {13, 0, 12};
Circle(14) = {12, 0, 11};

// Exterior boundary of the layer
Point(21) = {0, -(Rdom+Rpml), 0};
Point(22) = { (Rdom+Rpml), 0, 0};
Point(23) = {0,  (Rdom+Rpml), 0};
Point(24) = {-(Rdom+Rpml), 0, 0};
Circle(21) = {21, 0, 24};
Circle(22) = {24, 0, 23};
Circle(23) = {23, 0, 22};
Circle(24) = {22, 0, 21};

// Domain (inside disk)
Curve Loop(1) = {-1, -2, -3, -4};
Plane Surface(1) = {1};

// Domain (outside disk)
Curve Loop(2) = {1, 2, 3, 4, -11, -12, -13, -14};
Plane Surface(2) = {2};

// Layer
Curve Loop(3) = {11, 12, 13, 14, -21, -22, -23, -24};
Plane Surface(3) = {3};

// Definition of the physical regions
Physical Curve(201) = {21, 22, 23, 24}; // Exterior boundary of the layer
Physical Surface(301) = {1};            // Domain (inside disk)
Physical Surface(302) = {2};            // Domain (outside disk)
Physical Surface(303) = {3};            // Layer