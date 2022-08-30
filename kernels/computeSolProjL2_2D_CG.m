function [solP, sysP] = computeSolProjL2_2D_CG(mesh, dofm)

matP = sparse(dofm.numDofTRI, dofm.numDofTRI);
rhsP = zeros(dofm.numDofTRI, 1);

% Quadrature
degreeQ = 4*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);

% Shape functions
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
    
    % Elemental matrices
    matPel = shapeOrQ' * weights * shapeOrQ * detJdxdu;
    rhsPel = shapeOrQ' * weights * solQ * detJdxdu;
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matP(dof,dof) = matP(dof,dof) + matPel;
    rhsP(dof) = rhsP(dof) + rhsPel;
    
end

% Solution
solP = matP\rhsP;

sysP.matP = matP;
sysP.rhsP = rhsP;

end