function [matM, matK, matDX, matDY, matS1, matS2, matS3, matS4, dofS1, dofS2, dofS3, dofS4] = buildMatrixGlo2D_CG(mesh, dofm)

% -------------------------------------------------------------------------
% Compute volume matrices
% -------------------------------------------------------------------------

% Quadrature
degreeQ = 2*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);

% Shape functions (f, dfdu, dfdv)
funQ = functionsShapeTRI(uQ, vQ, dofm.degree);
[funDuQ, funDvQ] = functionsShapeDerTRI(uQ, vQ, dofm.degree);

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
    
    % Shape functions (dfdx, dfdy)
    funDxQ = funDuQ * Jdudx(1,1) + funDvQ * Jdudx(2,1);
    funDyQ = funDuQ * Jdudx(1,2) + funDvQ * Jdudx(2,2);
    
    % Elemental matrices
    matMel = zeros(dofm.numDofPerTRI,dofm.numDofPerTRI);
    matKel = zeros(dofm.numDofPerTRI,dofm.numDofPerTRI);
    matDXel = zeros(dofm.numDofPerTRI,dofm.numDofPerTRI);
    matDYel = zeros(dofm.numDofPerTRI,dofm.numDofPerTRI);
    for i=1:dofm.numDofPerTRI
        for j=1:dofm.numDofPerTRI
            matMel(i,j)  = weights' * (funQ(:,i) .* funQ(:,j)) * detJdxdu;
            matKel(i,j)  = weights' * (funDxQ(:,i) .* funDxQ(:,j) + funDyQ(:,i) .* funDyQ(:,j)) * detJdxdu;
            matDXel(i,j) = weights' * (funDxQ(:,i) .* funQ(:,j)) * detJdxdu;
            matDYel(i,j) = weights' * (funDyQ(:,i) .* funQ(:,j)) * detJdxdu;
        end
    end
    
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