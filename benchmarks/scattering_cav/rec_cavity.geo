//Mesh.MshFileVersion = 2.2;

Printf("Parameters : L = %f, L_PML = %f\n", L, L_PML);

// Points de la frontiere exterieure
Point(1) = {-(L+L_PML), -(L+L_PML), 0, 1.0};
Point(2) = {(L+L_PML), -(L+L_PML), 0, 1.0};
Point(3) = {(L+L_PML), (L+L_PML), 0, 1.0};
Point(4) = {-(L+L_PML), (L+L_PML), 0, 1.0};

// Points de la frontiere domaine-PML
Point(5) = {-L, -L, 0};
Point(6) = {L, -L, 0};
Point(7) = {L, L, 0};
Point(8) = {-L, L, 0};

// Frontiere exterieure
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Curve Loop(1) = {1, 2, 3, 4};

// Frontiere domaine-PML
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Curve Loop(2) = {7, 8, 5, 6};

// Points de la cavite rectangulaire
Point(9) = {0.75, -0.3, 0};
Point(10) = {-0.75, -0.3, 0};
Point(11) = {-0.75, -0.2, 0};
Point(12) = {0.55, -0.2, 0};
Point(13) = {0.55, 0.2, 0};
Point(14) = {-0.75, 0.2, 0};
Point(15) = {-0.75, 0.3, 0};
Point(16) = {0.75, 0.3, 0};

// Frontiere de la cavite rectangulaire
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 9};
Curve Loop(3) = {9, 10, 11, 12, 13, 14, 15, 16};

// Definition des surfaces du domaine et PML
Plane Surface(1) = {3, 2}; // domaine sans PML
Plane Surface(2) = {1, 2}; // PML

// Definition des elements physiques
Physical Curve("bd_pml", 1) = {1,2,3,4};
Physical Curve("obstacle", 2) = {9, 10, 11, 12, 13, 14, 15, 16};
Physical Surface(1) = {1};
Physical Surface("pml", 2) = {2};
