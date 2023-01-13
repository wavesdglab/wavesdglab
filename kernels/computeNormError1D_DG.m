function [normErr, normSol] = computeNormError1D_DG(mesh, dofm, mySol, vecSol)

Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeFunc = functionsShape1D(nodes,dofm.degree);

normSol2 = 0;
normErr2 = 0;
for e=1:mesh.numE
    
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2-coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    % Building the approximate solution
    numQ = zeros(1,Q);
    for n=1:dofm.numDofPerE
        numQ = numQ + vecSol(dofm.locToGlo(e,n)) * (shapeFunc(:,n))';
    end
    
    % Building the reference solution
    refQ = mySol(coordGlo);
    
    % Building the error
    errQ = numQ(:) - refQ(:);
    
    % Compute the errors
    normSol2  = normSol2  + weights(:)' * (refQ.*conj(refQ)) * (length/2) ;
    normErr2 = normErr2 + weights(:)' * (errQ.*conj(errQ)) * (length/2) ;
end

normSol = sqrt(normSol2);
normErr = sqrt(normErr2)/normSol;

end