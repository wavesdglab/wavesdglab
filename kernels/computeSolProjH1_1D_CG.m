function vecSolProjH1 = computeSolProjH1_1D_CG(mesh, dofm, mySol, mySolDer)

% Build matrix of the system

[matM, matK, matD] = buildMatrixGlo1D_CG(mesh, dofm);
matA = matK + matM;

% Set Dirichlet BC (left)

% matM(1,:) = 0;
% matM(1,1) = 1;

% Set Dirichlet BC (right)

% matM(mesh.numV,:) = 0;
% matM(mesh.numV,mesh.numV) = 1;

% Build RHS vector

Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeFunc = functionsShape1D(nodes,dofm.degree);
shapeFuncDer = functionsShapeDer1D(nodes,dofm.degree);

vecRhs = zeros(dofm.numDof,1);
for e=1:mesh.numE
    
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    solRef = mySol(coordGlo);
    derRef = mySolDer(coordGlo);
    %solRef = mySol(coordGlo) - mySolDerDer(coordGlo);
    
    glo = dofm.locToGlo(e,:);
    for i=1:dofm.numDofPerE
        vecRhs(glo(i)) = vecRhs(glo(i)) + weights' * (shapeFunc(:,i) .* solRef) * (length/2);
        vecRhs(glo(i)) = vecRhs(glo(i)) + weights' * (shapeFuncDer(:,i) .* derRef);
    end
end

% Set Neumann BC (left)

%vecRhs(1) = -mySolDer(0);

% Set Neumann BC (right)

%vecRhs(mesh.numV) = mySolDer(mesh.coordV(mesh.numV));

% Compute solution

vecSolProjH1 = matA\vecRhs;

end