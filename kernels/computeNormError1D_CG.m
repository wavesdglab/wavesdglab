function [errorL2sol, errorL2der, errorH1sol] = computeNormError1D_CG(mesh, dofm, mySol, mySolDer, vecSol)

errorL2sol2 = 0;
errorL2der2 = 0;

%Q = ceil((degree^2+1)/2);
Q = 16;
[nodes, weights] = quadratureGauss1D(Q);
shapeFunc = functionsShape1D(nodes,dofm.degree);
shapeFuncDer = functionsShapeDer1D(nodes,dofm.degree);

for e=1:mesh.numE
    
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2-coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    % Building the approximate solution (and its derivative)
    valNum = zeros(1,Q);
    valNumDer = zeros(1,Q);
    for n=1:dofm.numDofPerE
        valNum    = valNum    + vecSol(dofm.locToGlo(e,n)) * (shapeFunc(:,n)   )';
        valNumDer = valNumDer + vecSol(dofm.locToGlo(e,n)) * (shapeFuncDer(:,n))' * (2/length);
    end
    
    % Building the reference solution (and its derivative)
    valRef = mySol(coordGlo);
    valRefDer = mySolDer(coordGlo);
    
    % Building the error
    valErr = valNum(:) - valRef(:);
    valErrDer = valNumDer(:) - valRefDer(:);
    
    % Compute the errors
    errorL2sol2 = errorL2sol2 + weights(:)' * (valErr .* conj(valErr)) * (length/2) ;
    errorL2der2 = errorL2der2 + weights(:)' * (valErrDer .* conj(valErrDer)) * (length/2) ;
end

errorL2sol = sqrt(errorL2sol2);
errorL2der = sqrt(errorL2der2);
errorH1sol = sqrt(errorL2sol2 + errorL2der2);

end