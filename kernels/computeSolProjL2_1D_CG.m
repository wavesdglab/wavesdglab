function vecSolProjL2 = computeSolProjL2_1D_CG(mesh, dofm)

% Build matrix of the system
matElemM = buildMatrixElem1D(dofm.degree);
matMgg = sparse(dofm.numDofGam, dofm.numDofGam);
matMgi = sparse(dofm.numDofGam, dofm.numDofInt);
matMig = sparse(dofm.numDofInt, dofm.numDofGam);
matMii = sparse(dofm.numDofInt, dofm.numDofInt);
for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    
    matLocM = matElemM * length/2;
    
    glo = dofm.locToGlo(e,:);
    gloG = glo(1:2);
    gloI = glo(3:end)-dofm.numDofGam;
    
    matMgg(gloG,gloG) = matMgg(gloG,gloG) + matLocM(1:2,1:2);
    matMgi(gloG,gloI) = matMgi(gloG,gloI) + matLocM(1:2,3:end);
    matMig(gloI,gloG) = matMig(gloI,gloG) + matLocM(3:end,1:2);
    matMii(gloI,gloI) = matMii(gloI,gloI) + matLocM(3:end,3:end);
end
matM = [ matMgg matMgi ; matMig matMii ];

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

vecRhs = zeros(dofm.numDof,1);
for e=1:mesh.numE
    
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    refQ = mySol1D(coordGlo);
    
    glo = dofm.locToGlo(e,:);
    for i=1:dofm.numDofPerE
        vecRhs(glo(i)) = vecRhs(glo(i)) + weights' * (shapeFunc(:,i) .* refQ) * (length/2);
    end
end

% Compute solution

vecSolProjL2 = matM\vecRhs;

end