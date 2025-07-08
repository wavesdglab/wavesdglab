function solP = computeSolProjL2_1D_DG(mesh, dofm)

% Quadrature and shape functions
N = dofm.degree+1;
Q = ceil((2*(N-1)+1)/2);
Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeQ = functionsShapeLIN(nodes,dofm.degree);

% Compute (element-wise) L2-projection of the solution
solP = zeros(dofm.numDof,1);
for e=1:mesh.numE
    
    % Mapping
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    % Reference solution
    refQ = mySol(coordGlo);
    
    % Local matrix and RHS vector
    matLocP = shapeQ' * (weights .* shapeQ) * length/2;
    rhsLocP = shapeQ' * (weights .* refQ) * length/2;
    
    % Compute
    glo = dofm.locToGlo(e,:);
    solP(glo) = matLocP\rhsLocP;
end

end
