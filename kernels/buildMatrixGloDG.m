function [matM, matK, matDX, matDY] = buildMatrixGloDG(mesh, dofm)

matM = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Mass matrix
matK = sparse(dofm.numDofTRI,dofm.numDofTRI);   % Stiffness matrix
matDX = sparse(dofm.numDofTRI,dofm.numDofTRI);  % Differentiation matrix (x)
matDY = sparse(dofm.numDofTRI,dofm.numDofTRI);  % Differentiation matrix (y)

for tri=1:mesh.numTri
    
    % Elemental matrices
    verTri = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(verTri(1),:);
    V2 = mesh.coord(verTri(2),:);
    V3 = mesh.coord(verTri(3),:);
    [matMel, matKel, matDXel, matDYel] = buildMatrixElemTRI(V1,V2,V3,dofm.degree);
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matM(dof,dof) = matMel;
    matK(dof,dof) = matKel;
    matDX(dof,dof) = matDXel;
    matDY(dof,dof) = matDYel;
    
end

end