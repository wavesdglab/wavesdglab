// Dimensions
Lx = 67437.5;
Ly = 11943.75;
X0 = 0;
Y0 = -Ly;

// Parameters
X_SOU = Lx/2;
Y_SOU = -10;
X0 = 2*Lx/5;
Lx = Lx/5;
Ly = Ly/3;
Y0 = -Ly;

// ====================================================================================================
// Build DOMAIN
// ====================================================================================================

p1 = newp; Point(p1) = {X0 + Lx*0, Y0 + Ly*0, 0};
p2 = newp; Point(p2) = {X0 + Lx*1, Y0 + Ly*0, 0};
p3 = newp; Point(p3) = {X0 + Lx*1, Y0 + Ly*1, 0};
p4 = newp; Point(p4) = {X0 + Lx*0, Y0 + Ly*1, 0};

l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};

ll = newll; Line loop(ll) = {l1, l2, l3, l4};
s = news; Plane Surface(s) = {ll};

// ====================================================================================================
// Build SOURCE
// ====================================================================================================

p0 = newp; Point(p0) = {X_SOU, Y_SOU, 0};
Point {p0} In Surface {s};

// ====================================================================================================
// Generate MESH
// ====================================================================================================

// Physical/mesh parameters
// WAVENUMBER = 2*Pi*FREQ/C[X];
// LAMBDA = 2*Pi/WAVENUMBER;
// LC = LAMBDA/N_LAMBDA;
// LC = cMin/FREQ/N_LAMBDA;

Field[1] = Structured;
Field[1].FileName = "benchmarks/geophysics_BPmodel/_velocityGmsh.txt";
Field[1].TextFormat = 1;
Field[2] = MathEval;
Field[2].F = Sprintf("F1 * %g", 1/FREQ/N_LAMBDA);
Background Field = 2;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// ====================================================================================================
// Build PHYSICALS
// ====================================================================================================

TAG_SOU = 1000;
TAG_BND = 2000;
TAG_DOM = 3000;

Physical Point(TAG_SOU) = p0;
Physical Line(TAG_BND) = {l1, l2, l3, l4};
Physical Surface(TAG_DOM) = s;
