function [matM, matK, matDX, matDY, matS1, matS2, matS3, matS4, dofS1, dofS2, dofS3, dofS4] = buildMatrixGlo2D_CG(mesh, dofm)

% Quadrature
degreeQ = 2*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);

% Shape functions (f, dfdu, dfdv)
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);
[shapeDuQ, shapeDvQ] = functionsShapeDerTRI(uQ, vQ, dofm.degree);

% Global matrices
matM = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Mass matrix
matK = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Stiffness matrix
matDX = sparse(dofm.numDofTRI,dofm.numDofTRI);  % Differentiation matrix (x)
matDY = sparse(dofm.numDofTRI,dofm.numDofTRI);  % Differentiation matrix (y)

for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
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
    
    % Elemental matrices
    matMel  = shapeOrQ' * weights * shapeOrQ * detJdxdu;
    matKel  = (shapeDxQ' * weights * shapeDxQ + shapeDyQ' * weights * shapeDyQ ) * detJdxdu;
    matDXel = shapeDxQ' * weights * shapeOrQ * detJdxdu;
    matDYel = shapeDyQ' * weights * shapeOrQ * detJdxdu;
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matM(dof,dof) = matM(dof,dof) + matMel;
    matK(dof,dof) = matK(dof,dof) + matKel;
    matDX(dof,dof) = matDX(dof,dof) + matDXel;
    matDY(dof,dof) = matDY(dof,dof) + matDYel;
    
end

% -------------------------------------------------------------------------
% Compute boundary matrices
% -------------------------------------------------------------------------

matS1 = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Surfacic mass matrix on Gamma_1
matS2 = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Surfacic mass matrix on Gamma_2
matS3 = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Surfacic mass matrix on Gamma_3
matS4 = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Surfacic mass matrix on Gamma_4
dofS1 = zeros(dofm.numDofTRI,1);
dofS2 = zeros(dofm.numDofTRI,1);
dofS3 = zeros(dofm.numDofTRI,1);
dofS4 = zeros(dofm.numDofTRI,1);

for edgBnd=1:mesh.numEdgBnd
    edg = mesh.listEdgBnd(edgBnd);
    ver = mesh.mapEdgToVer(edg,:);
    
    % Elemental matrix
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    [Mel, ~, ~] = buildMatrixElemLIN(V1,V2,dofm.degree);
    
    % Matrix assembling
    switch mesh.tagEdgBnd(edgBnd)
        case 1
            matS1(ver,ver) = matS1(ver,ver) + Mel;
            dofS1(ver,1) = 1;
        case 2
            matS2(ver,ver) = matS2(ver,ver) + Mel;
            dofS2(ver,1) = 1;
        case 3
            matS3(ver,ver) = matS3(ver,ver) + Mel;
            dofS3(ver,1) = 1;
        case 4
            matS4(ver,ver) = matS4(ver,ver) + Mel;
            dofS4(ver,1) = 1;
    end
end

dofS1 = find(dofS1);
dofS2 = find(dofS2);
dofS3 = find(dofS3);
dofS4 = find(dofS4);

end