function vecSolProjL2 = computeSolProjL2_1D_DG(mesh, dofm)

% Quadrature and shape functions
Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeFunc = functionsShape1D(nodes,dofm.degree);

% Compute L2-projection of the solution

matElemM = buildMatrixElem1D(dofm.degree);
vecSolProjL2 = zeros(dofm.numDof,1);

for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    
    matLocM = matElemM * length/2;
    
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    solRef = mySol1D(coordGlo);
    
    glo = dofm.locToGlo(e,:);
    vecSolProjL2(glo) = matLocM \ shapeFunc' * (weights .* solRef) * (length/2);
end

end
