function [matM, matDX, matDY, matK] = buildMatrixElem2D(degree)

% Quadrature
degreeQ = 2*degree;
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);
weightsTriQ = sparse(1:size(weightsTriQ,1), 1:size(weightsTriQ,1), weightsTriQ);

% Shape functions
shapeTriQ = functionsShapeTRI(uTriQ, vTriQ, degree);
[shapeTriDuQ, shapeTriDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, degree);

V1 = [-1 -1];
V2 = [ 1 -1];
V3 = [-1  1];

Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
Jdudx = inv(Jdxdu);                 % [ du/dx du/dy ; dv/dx dv/dy ]
detJdxdu = abs(det(Jdxdu));

% Shape functions (f, dfdx, dfdy) with orientation
shapeDxQ = (shapeTriDuQ * Jdudx(1,1) + shapeTriDvQ * Jdudx(2,1));
shapeDyQ = (shapeTriDuQ * Jdudx(1,2) + shapeTriDvQ * Jdudx(2,2));

% Elemental matrices
matM  = shapeTriQ' * weightsTriQ * shapeTriQ * detJdxdu;
matDX = shapeDxQ'  * weightsTriQ * shapeTriQ * detJdxdu;
matDY = shapeDyQ'  * weightsTriQ * shapeTriQ * detJdxdu;
matK  = shapeDxQ'  * weightsTriQ * shapeDxQ * detJdxdu ...
      + shapeDyQ'  * weightsTriQ * shapeDyQ * detJdxdu;

end