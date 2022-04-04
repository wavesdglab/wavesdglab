function [matM, matK, matDX, matDY] = buildMatrixGlo2D_DG(mesh, dofm)

% Quadrature 2D
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
    matM(dof,dof) = matMel;
    matK(dof,dof) = matKel;
    matDX(dof,dof) = matDXel;
    matDY(dof,dof) = matDYel;
    
end

end