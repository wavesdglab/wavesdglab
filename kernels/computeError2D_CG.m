function [errorL2, errorH1, normL2, normH1] = computeError2D_CG(mesh, dofm, vecSol, vecRef)

degreeQ = 2*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
shapeFunc = functionsShapeTRI(uQ, vQ, dofm.degree);
[shapeFuncDu, shapeFuncDv] = functionsShapeDerTRI(uQ, vQ, dofm.degree);

normL2sol2 = 0;
normL2der2 = 0;
errorL2sol2 = 0;
errorL2der2 = 0;

for tri=1:mesh.numTri
    
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uQ, vQ, V1, V2, V3);
    Jdxdu = [(V3-V2)' (V1-V2)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    Jdudx = inv(Jdxdu);                 % [ du/dx du/dy ; dv/dx dv/dy ]
    detJdxdu = abs(det(Jdxdu));
    
    % Building the approximate solution (and derivatives)
    solQ   = shapeFunc * vecSol(dofm.locToGloTRI(tri,:));
    solDxQ = (shapeFuncDu * Jdudx(1,1) + shapeFuncDv * Jdudx(2,1)) * vecSol(dofm.locToGloTRI(tri,:));
    solDyQ = (shapeFuncDu * Jdudx(1,2) + shapeFuncDv * Jdudx(2,2)) * vecSol(dofm.locToGloTRI(tri,:));
    
    % Building the reference solution (and derivatives) for a reference vector
    if (exist('vecRef','var'))
        refQ   = shapeFunc * vecRef(dofm.locToGloTRI(tri,:));
        refDxQ = (shapeFuncDu * Jdudx(1,1) + shapeFuncDv * Jdudx(2,1)) * vecRef(dofm.locToGloTRI(tri,:));
        refDyQ = (shapeFuncDu * Jdudx(1,2) + shapeFuncDv * Jdudx(2,2)) * vecRef(dofm.locToGloTRI(tri,:));
    else
        [refQ, refDxQ, refDyQ] = mySol(xQ, yQ);
    end
    
    % Building the error
    errQ   = solQ(:)   - refQ(:);
    errDxQ = solDxQ(:) - refDxQ(:);
    errDyQ = solDyQ(:) - refDyQ(:);
    
    % Compute the errors
    normL2sol2 = normL2sol2 + weights(:)' * (refQ .* conj(refQ)) * detJdxdu;
    normL2der2 = normL2der2 + weights(:)' * (refDxQ .* conj(refDxQ) + refDyQ .* conj(refDyQ)) * detJdxdu;
    errorL2sol2 = errorL2sol2 + weights(:)' * (errQ .* conj(errQ)) * detJdxdu;
    errorL2der2 = errorL2der2 + weights(:)' * (errDxQ .* conj(errDxQ) + errDyQ .* conj(errDyQ)) * detJdxdu;
end

normL2 = sqrt(normL2sol2);
normH1 = sqrt(normL2sol2 + normL2der2);
errorL2 = sqrt(errorL2sol2)/normL2;
errorH1 = sqrt(errorL2sol2 + errorL2der2)/normH1;

end