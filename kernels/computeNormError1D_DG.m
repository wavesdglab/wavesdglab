function errorL2sol = computeNormError1D_DG(mesh, dofm, mySol, vecSol)

Q = 16;
[nodes, weights] = quadratureGauss1D(Q);
shapeFunc = functionsShape1D(nodes,dofm.degree);

errorL2sol2 = 0;
for e=1:mesh.numE
    
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2-coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    % Building the approximate solution
    solNum = zeros(1,Q);
    for n=1:dofm.numDofPerE
        solNum = solNum + vecSol(dofm.locToGlo(e,n)) * (shapeFunc(:,n))';
    end
    
    % Building the reference solution
    solRef = mySol(coordGlo);
    
    % Building the error
    errNum = solNum(:) - solRef(:);
    
    % Compute the errors
    errorL2sol2 = errorL2sol2 + weights(:)' * (errNum .* conj(errNum)) * (length/2) ;
end

errorL2sol = sqrt(errorL2sol2);

end