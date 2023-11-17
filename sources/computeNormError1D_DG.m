% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [errorL2, normL2] = computeNormError1D_DG(mesh, dofm, vecSol)

% Quadrature and shape functions
Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeQ = functionsShape1D(nodes,dofm.degree);

normRef2 = 0;
normErr2 = 0;
for e=1:mesh.numE
    
    % Mapping
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2-coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    % Numerical solution, reference solution and error
    solQ = shapeQ * vecSol(dofm.locToGlo(e,:));
%     refQ = mySol1D(coordGlo);
    refQ = mySol1D_heterogeneous(coordGlo);
    errQ = solQ(:) - refQ(:);
    
    % Norms
    normRef2 = normRef2 + weights' * (refQ.*conj(refQ)) * (length/2) ;
    normErr2 = normErr2 + weights' * (errQ.*conj(errQ)) * (length/2) ;
end

normL2  = sqrt(normRef2);
errorL2 = sqrt(normErr2)/normL2;

end