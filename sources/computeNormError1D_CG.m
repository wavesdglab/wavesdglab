% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [errorL2, errorH1, normL2, normH1] = computeNormError1D_CG(mesh, dofm, vecSol)

% Quadrature and shape functions
Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeQ = functionsShapeLIN(nodes,dofm.degree);
shapeDerQ = functionsShapeDerLIN(nodes,dofm.degree);

normL2sol2 = 0;
normL2der2 = 0;
errorL2sol2 = 0;
errorL2der2 = 0;
for e=1:mesh.numE
    
    % Mapping
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2-coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    % Numerical solution, reference solution, error + all the derivatives
    solQ = shapeQ * vecSol(dofm.locToGlo(e,:));
    solDerQ = shapeDerQ * vecSol(dofm.locToGlo(e,:)) * (2/length);
    [refQ, refDerQ] = mySol(coordGlo);
    errQ = solQ - refQ;
    errDerQ = solDerQ - refDerQ;
    
    % Norms
    normL2sol2  = normL2sol2  + weights' * (refQ .* conj(refQ)) * (length/2);
    normL2der2  = normL2der2  + weights' * (refDerQ .* conj(refDerQ)) * (length/2);
    errorL2sol2 = errorL2sol2 + weights' * (errQ .* conj(errQ)) * (length/2) ;
    errorL2der2 = errorL2der2 + weights' * (errDerQ .* conj(errDerQ)) * (length/2) ;
end

normL2 = sqrt(normL2sol2);
normH1 = sqrt(normL2sol2 + normL2der2);
errorL2 = sqrt(errorL2sol2)/normL2;
errorH1 = sqrt(errorL2sol2 + errorL2der2)/normH1;

end