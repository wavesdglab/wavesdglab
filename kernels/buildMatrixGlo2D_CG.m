function [matM, matK, matDX, matDY, matS1, matS2, matS3, matS4, dofS1, dofS2, dofS3, dofS4] = buildMatrixGlo2D_CG(mesh, dofm)

% -------------------------------------------------------------------------
% Compute volume matrices
% -------------------------------------------------------------------------

matM = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Mass matrix
matK = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Stiffness matrix
matDX = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Differentiation matrix (x)
matDY = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Differentiation matrix (y)

for tri=1:mesh.numTri
    
    % Elemental matrices
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [Mel, Kel, DXel, DYel] = buildMatrixElemTRI(V1,V2,V3,dofm.degree);
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matK(dof,dof) = matK(dof,dof) + Kel;
    matM(dof,dof) = matM(dof,dof) + Mel;
    matDX(dof,dof) = matDX(dof,dof) + DXel;
    matDY(dof,dof) = matDY(dof,dof) + DYel;
    
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