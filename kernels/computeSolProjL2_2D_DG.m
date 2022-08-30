function [solP, sysP] = computeSolProjL2_2D_DG(mesh, dofm)

matPx = zeros(mesh.numTri*dofm.numDofPerTRI, dofm.numDofPerTRI);
matPy = zeros(mesh.numTri*dofm.numDofPerTRI, dofm.numDofPerTRI);
matPv = zeros(mesh.numTri*dofm.numDofPerTRI, dofm.numDofPerTRI);
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
    matPx(dof,:) = dof'*ones(1,dofm.numDofPerTRI);
    matPy(dof,:) = ones(dofm.numDofPerTRI,1)*dof;
    matPv(dof,:) = matPel;
    rhsP(dof) = rhsPel;
    
end

matP = sparse(matPx, matPy, matPv, dofm.numDofTRI, dofm.numDofTRI);

% Solution
solP = matP\rhsP;

sysP.matP = matP;
sysP.rhsP = rhsP;

end