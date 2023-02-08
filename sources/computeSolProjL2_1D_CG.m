% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function solP = computeSolProjL2_1D_CG(mesh, dofm)

% Quadrature and shape functions
N = dofm.degree+1;
Q = ceil((2*(N-1)+1)/2);
Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeQ = functionsShape1D(nodes,dofm.degree);

% Build matrix and RHS vector
matP = sparse(dofm.numDof, dofm.numDof);
rhsP = zeros(dofm.numDof,1);
for e=1:mesh.numE
    
    % Mapping
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    % Reference solution
    refQ = mySol1D(coordGlo);
    
    % Local matrix and RHS vector
    matLocP = shapeQ' * (weights .* shapeQ) * length/2;
    rhsLocP = shapeQ' * (weights .* refQ) * length/2;
    
    % Assembling
    glo = dofm.locToGlo(e,:);
    matP(glo,glo) = matP(glo,glo) + matLocP;
    rhsP(glo) = rhsP(glo) + rhsLocP;
end

% Set Dirichlet BC (left)
% matP(1,:) = 0;
% matP(1,1) = 1;

% Set Dirichlet BC (right)
% matP(mesh.numV,:) = 0;
% matP(mesh.numV,mesh.numV) = 1;

% Compute L2-projection of the solution
solP = matP\rhsP;

end