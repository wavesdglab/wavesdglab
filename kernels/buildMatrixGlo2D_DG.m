function [matM, matK, matDX, matDY, rhsP] = buildMatrixGlo2D_DG(mesh, dofm)

% Quadrature
degreeQ = 2*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);

% Shape functions (f, dfdu, dfdv)
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);
[shapeDuQ, shapeDvQ] = functionsShapeDerTRI(uQ, vQ, dofm.degree);

% Global matrices
matM = sparse(dofm.numDofTRI, dofm.numDofTRI);   % Mass matrix
matK = sparse(dofm.numDofTRI, dofm.numDofTRI);   % Stiffness matrix
matDX = sparse(dofm.numDofTRI, dofm.numDofTRI);  % Differentiation matrix (x)
matDY = sparse(dofm.numDofTRI, dofm.numDofTRI);  % Differentiation matrix (y)
rhsP = zeros(dofm.numDofTRI, 1);

for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uQ, vQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    Jdudx = inv(Jdxdu);                 % [ du/dx du/dy ; dv/dx dv/dy ]
    detJdxdu = abs(det(Jdxdu));
    
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
    
    % Shape functions (f, dfdx, dfdy) with orientation
    shapeOrQ = shapeQ * orientation;
    shapeDxQ = (shapeDuQ * Jdudx(1,1) + shapeDvQ * Jdudx(2,1)) * orientation;
    shapeDyQ = (shapeDuQ * Jdudx(1,2) + shapeDvQ * Jdudx(2,2)) * orientation;
    
    % RHS function
    [~, ~, ~, rhsQ] = mySol(xQ, yQ);
    
    % Elemental matrices
    matMel = shapeOrQ' * weights * shapeOrQ * detJdxdu;
    matKel = shapeDxQ' * weights * shapeDxQ * detJdxdu + shapeDyQ' * weights * shapeDyQ * detJdxdu;
    matDXel = shapeDxQ' * weights * shapeOrQ * detJdxdu;
    matDYel = shapeDyQ' * weights * shapeOrQ * detJdxdu;
    vecRHSel = shapeOrQ' * weights * rhsQ * detJdxdu;
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matM(dof,dof) = matMel;
    matK(dof,dof) = matKel;
    matDX(dof,dof) = matDXel;
    matDY(dof,dof) = matDYel;
    rhsP(dof) = vecRHSel;
    
end

end