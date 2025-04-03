SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
Circle(2) = {-0, -0, 0, 1, 0, 2*Pi};
Curve Loop(1) = {1};
Plane Surface(1) = {1};
Physical Curve(3) = {1};
Physical Surface(4) = {1};

X_SOU = -M;
Y_SOU = 0;
p0 = newp; Point(p0) = {X_SOU, Y_SOU, 0};
Point {p0} In Surface {1};

TAG_SOU = 1;
Physical Point(TAG_SOU) = p0;

Point(100) = {X_SOU, Y_SOU, 0, 0.03};

\\Field[1] = Distance;
\\Field[1].NodesList = {100};
\\Field[2] = Threshold;
\\Field[2].IField = 1;
\\Field[2].LcMin = hmin;
\\Field[2].LcMax = h;
\\Field[2].DistMin = 0.05;
\\Field[2].DistMax = 0.05;
\\Background Field = 2;