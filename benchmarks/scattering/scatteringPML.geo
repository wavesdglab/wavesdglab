//Mesh.MshFileVersion = 2.2;

Printf("Parameters : L = %f, L_PML = %f, R_disk = %f\n", L, L_PML, R_disk);

// Points de la frontiere exterieure
Point(1) = {-(L+L_PML), -(L+L_PML), 0, 1.0};
Point(2) = {(L+L_PML), -(L+L_PML), 0, 1.0};
Point(3) = {(L+L_PML), (L+L_PML), 0, 1.0};
Point(4) = {-(L+L_PML), (L+L_PML), 0, 1.0};

// Points de la frontiere du disque
Point(5) = {0, -R_disk, 0};
Point(6) = {R_disk, 0, 0};
Point(7) = {0, R_disk, 0};
Point(8) = {-R_disk, 0, 0};

// Centre du disque
Point(9) = {0, 0, 0};

// Points de la frontiere domaine-PML
Point(10) = {-L, -L, 0};
Point(11) = {L, -L, 0};
Point(12) = {L, L, 0};
Point(13) = {-L, L, 0};

// Frontiere exterieure
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Curve Loop(1) = {1, 2, 3, 4};

// Disque
Circle(5) = {6, 9, 7};
Circle(6) = {7, 9, 8};
Circle(7) = {8, 9, 5};
Circle(8) = {5, 9, 6};
Curve Loop(2) = {5, 6, 7, 8};

// Frontiere domaine-PML
Line(9) = {10, 11};
Line(10) = {11, 12};
Line(11) = {12, 13};
Line(12) = {13, 10};
Curve Loop(3) = {11, 12, 9, 10};

// Definition des surfaces du domaine et PML
Plane Surface(1) = {3, 2};
Plane Surface(2) = {1, 3};

// Definition des elements physiques
Physical Curve("bd_pml", 1) = {1,2,3,4};
Physical Curve("bd_obstacle", 2) = {5,6,7,8};
Physical Surface("domain", 1) = {1};
Physical Surface("pml", 2) = {2};
