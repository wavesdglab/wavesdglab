function [solP, matP, rhsP] = computeSolProj2D_CG(mesh, dofm)

matP = sparse(dofm.numDofTRI, dofm.numDofTRI);
rhsP = zeros(dofm.numDofTRI, 1);

% Quadrature
degreeQ = 4*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);

% Shape functions (f, dfdu, dfdv)
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);

for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uQ, vQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    detJdxdu = abs(det(Jdxdu));
    
    % RHS function
    solQ = mySol(xQ, yQ);
    
    % Elemental matrices
    matPel = shapeQ' * weights * shapeQ * detJdxdu;
    rhsPel = shapeQ' * weights * solQ * detJdxdu;
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matP(dof,dof) = matP(dof,dof) + matPel;
    rhsP(dof) = rhsP(dof) + rhsPel;
    
end

% Solution
solP = matP\rhsP;

end