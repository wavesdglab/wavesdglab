Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {0, 1, 0};
Point(4) = {-1, 0, 0};
Point(5) = {0, -1, 0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Physical Curve(3) = {1,2,3,4};
Physical Surface(4) = {1};

p0 = newp; Point(p0) = {xSou, ySou, 0};
Point {p0} In Surface {1};

TAG_SOU = 1;
Physical Point(TAG_SOU) = p0;

Field[1] = Distance;
Field[1].NodesList = {p0};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = h/5;
Field[2].LcMax = h*2;
Field[2].DistMin = h;
Field[2].DistMax = h*5;
Field[3] = MathEval;
Field[3].F = Sprintf("(-((x-(%g))/max(sqrt((x-(%g))^2+y^2),0.01))*(%g) + sqrt(1 + (y/max(sqrt((x-(%g))^2+y^2),0.01))^2*(%g)^2))*(%g)",xSou,xSou,xSou,xSou,xSou,h);
Field[3].F = Sprintf("(-((x-(%g))/max(sqrt((x-(%g))^2+y^2),0.01))*(%g) + sqrt(1 + (y/max(sqrt((x-(%g))^2+y^2),0.01))^2*(%g)^2))*(%g)",xSou,xSou,xSou,xSou,xSou,h);
Field[4] = Min;
Field[4].FieldsList = {2,3};
Background Field = 4;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 5;

//x = -1:0.02:1;
//y = -1:0.02:1;
//[X,Y] = meshgrid(x,y);
//x0 = -0.25;
//D = sqrt((X-x0).^2+Y.^2);
//surf(X,Y,D)
//SIN = Y./D;
//COS = (X-x0)./D;
//DIS = -COS*x0 + sqrt(1+SIN.^2*x0^2);
//surf(X,Y,DIS)