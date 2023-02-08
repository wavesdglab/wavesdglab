% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [sol] = computeSolProjL2_2D_DG(mesh, dofm)

% Quadrature and shape functions
degreeQ = 4*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);

% Build matrix and RHS vector
matPx = zeros(mesh.numTri*dofm.numDofPerTRI, dofm.numDofPerTRI);
matPy = zeros(mesh.numTri*dofm.numDofPerTRI, dofm.numDofPerTRI);
matPv = zeros(mesh.numTri*dofm.numDofPerTRI, dofm.numDofPerTRI);
rhsU  = zeros(dofm.numDofTRI,1);
rhsVx = zeros(dofm.numDofTRI,1);
rhsVy = zeros(dofm.numDofTRI,1);
for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uQ, vQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    detJdxdu = abs(det(Jdxdu));
    
    % Reference solution
    [refQ, ~, ~, ~, refVxQ, refVyQ] = mySol(xQ, yQ);
    
    % Orientation
    orientation = ones(dofm.numDofPerTRI,1);
    if(ver(1) > ver(2))
        orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(ver(2) > ver(3))
        orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(ver(3) > ver(1))
        orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    orientation = sparse(1:dofm.numDofPerTRI, 1:dofm.numDofPerTRI, orientation);
    
    % Shape functions with orientation
    shapeOrQ = shapeQ * orientation;
    
    % Local matrix and RHS vector
    matPel = shapeOrQ' * weights * shapeOrQ * detJdxdu;
    rhsUel = shapeOrQ' * weights * refQ * detJdxdu;
    rhsVxel = shapeOrQ' * weights * refVxQ * detJdxdu;
    rhsVyel = shapeOrQ' * weights * refVyQ * detJdxdu;
    
    % Assembling
    dof = dofm.locToGloTRI(tri,:);
    matPx(dof,:) = dof'*ones(1,dofm.numDofPerTRI);
    matPy(dof,:) = ones(dofm.numDofPerTRI,1)*dof;
    matPv(dof,:) = matPel;
    rhsU(dof)  = rhsUel;
    rhsVx(dof) = rhsVxel;
    rhsVy(dof) = rhsVyel;
end

matP = sparse(matPx, matPy, matPv, dofm.numDofTRI, dofm.numDofTRI);

% Solution
solU  = matP\rhsU;
solVx = matP\rhsVx;
solVy = matP\rhsVy;

sol = [solU ; solVx ; solVy];

end