function vecSolProjL2 = computeSolProjL2_1D_DG(mesh, dofm, mySol)

% Build matrix of the system

[matM, matK, matD] = buildMatrixGlo1D_DG(mesh, dofm);

% Build RHS vector

Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeFunc = functionsShape1D(nodes,dofm.degree);

vecRhs = zeros(dofm.numDof,1);
for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    solRef = mySol(coordGlo);
    
    glo = dofm.locToGlo(e,:);
    vecRhs(glo) = vecRhs(glo) + shapeFunc' * (weights .* solRef) * (length/2);
end

% Compute solution

vecSolProjL2 = matM\vecRhs;

end
