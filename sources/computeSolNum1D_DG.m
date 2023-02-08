% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solA, sysA] = computeSolNum1D_DG(mesh, dofm, theta, tau, PREC)

global k BCLeft BCRight

% -------------------------------------------------------------------------
% Quadrature and shape functions
% -------------------------------------------------------------------------

Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeQ = functionsShape1D(nodes,dofm.degree);
shapeDerQ = functionsShapeDer1D(nodes,dofm.degree);
matMelem = shapeQ' * (weights .* shapeQ);
matDelem = shapeQ' * (weights .* shapeDerQ);

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

matM = sparse(dofm.numDof, dofm.numDof);
matD = sparse(dofm.numDof, dofm.numDof);
rhsP = zeros(dofm.numDof,1);
rhsU = zeros(dofm.numDof,1);

for e=1:mesh.numE
    
    % Mapping
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    % Local RHS vector
    [~, ~, ~, ~, souP, souU] = mySol1D(coordGlo);
    rhsPloc = (shapeQ .* souP).' * weights * (length/2);
    rhsUloc = (shapeQ .* souU).' * weights * (length/2);
    
    % Local matrices
    matMloc = matMelem * length/2;
    matDloc = matDelem;
    
    % Assembling
    glo = dofm.locToGlo(e,:);
    matM(glo,glo) = matM(glo,glo) + matMloc;
    matD(glo,glo) = matD(glo,glo) + matDloc;
    rhsP(glo) = rhsP(glo) + rhsPloc;
    rhsU(glo) = rhsU(glo) + rhsUloc;
end

matA = [ -1i*k*matM -matD' ; -matD' -1i*k*matM ];
rhsA = [ rhsP ; rhsU ];

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

[solPL, ~, solUL] = mySol1D(0);
[solPR, ~, solUR] = mySol1D(mesh.coordV(mesh.numV));

for e=1:mesh.numE
    
    % Left
    idIntP = dofm.locToGlo(e,1);
    idIntU = dofm.locToGlo(e,1) + dofm.numDof;
    
    if(e > 1)
        idExtP = dofm.locToGlo(e-1,2);
        idExtU = dofm.locToGlo(e-1,2) + dofm.numDof;
        matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau * theta;
        matA(idIntP,idIntU) = matA(idIntP,idIntU) - 0.5;
        matA(idIntP,idExtP) = matA(idIntP,idExtP) - 0.5*tau * theta;
        matA(idIntP,idExtU) = matA(idIntP,idExtU) - 0.5;
        matA(idIntU,idIntP) = matA(idIntU,idIntP) - 0.5;
        matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau * theta;
        matA(idIntU,idExtP) = matA(idIntU,idExtP) - 0.5;
        matA(idIntU,idExtU) = matA(idIntU,idExtU) - 0.5/tau * theta;
    else
        switch BCLeft
            case 'PER'
                idExtP = dofm.locToGlo(mesh.numE,2);
                idExtU = dofm.locToGlo(mesh.numE,2) + dofm.numDof;
                matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau * theta;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) - 0.5;
                matA(idIntP,idExtP) = matA(idIntP,idExtP) - 0.5*tau * theta;
                matA(idIntP,idExtU) = matA(idIntP,idExtU) - 0.5;
                matA(idIntU,idIntP) = matA(idIntU,idIntP) - 0.5;
                matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau * theta;
                matA(idIntU,idExtP) = matA(idIntU,idExtP) - 0.5;
                matA(idIntU,idExtU) = matA(idIntU,idExtU) - 0.5/tau * theta;
            case 'DIR'
                matA(idIntP,idIntP) = matA(idIntP,idIntP) + tau * theta;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) - 1;
                rhsA(idIntP) = rhsA(idIntP) + solPL*tau * theta;
                rhsA(idIntU) = rhsA(idIntU) + solPL;
            case 'NEU'
                matA(idIntU,idIntP) = matA(idIntU,idIntP) - 1;
                matA(idIntU,idIntU) = matA(idIntU,idIntU) + 1/tau * theta;
                rhsA(idIntP) = rhsA(idIntP) + solUL;
                rhsA(idIntU) = rhsA(idIntU) + solUL/tau * theta;
            case 'ABC'
                matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) - 0.5;
                matA(idIntU,idIntP) = matA(idIntU,idIntP) - 0.5;
                matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5;
            otherwise
                warning('Error - No valid BC has been set on the left.')
        end
    end
    
    % Right
    idIntP = dofm.locToGlo(e,2);
    idIntU = dofm.locToGlo(e,2) + dofm.numDof;
    
    if(e < mesh.numE)
        idExtP = dofm.locToGlo(e+1,1);
        idExtU = dofm.locToGlo(e+1,1) + dofm.numDof;
        matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau * theta;
        matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5;
        matA(idIntP,idExtP) = matA(idIntP,idExtP) - 0.5*tau * theta;
        matA(idIntP,idExtU) = matA(idIntP,idExtU) + 0.5;
        matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5;
        matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau * theta;
        matA(idIntU,idExtP) = matA(idIntU,idExtP) + 0.5;
        matA(idIntU,idExtU) = matA(idIntU,idExtU) - 0.5/tau * theta;
    else
        switch BCRight
            case 'PER'
                idExtP = dofm.locToGlo(1,1);
                idExtU = dofm.locToGlo(1,1) + dofm.numDof;
                matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau * theta;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5;
                matA(idIntP,idExtP) = matA(idIntP,idExtP) - 0.5*tau * theta;
                matA(idIntP,idExtU) = matA(idIntP,idExtU) + 0.5;
                matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5;
                matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau * theta;
                matA(idIntU,idExtP) = matA(idIntU,idExtP) + 0.5;
                matA(idIntU,idExtU) = matA(idIntU,idExtU) - 0.5/tau * theta;
            case 'DIR'
                matA(idIntP,idIntP) = matA(idIntP,idIntP) + tau * theta;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) + 1;
                rhsA(idIntP) = rhsA(idIntP) + solPR*tau * theta;
                rhsA(idIntU) = rhsA(idIntU) - solPR;
            case 'NEU'
                matA(idIntU,idIntP) = matA(idIntU,idIntP) + 1;
                matA(idIntU,idIntU) = matA(idIntU,idIntU) + 1/tau * theta;
                rhsA(idIntP) = rhsA(idIntP) - solUR;
                rhsA(idIntU) = rhsA(idIntU) + solUR/tau * theta;
            case 'ABC'
                matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5;
                matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5;
                matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5;
            otherwise
                warning('Error - No valid BC has been set on the left.')
        end
    end
    
end

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

% Preconditioning
switch PREC
    case 'PrecMass'
        matP = [matM zeros(size(matM,1)) ; zeros(size(matM,1)) matM];
    case 'PrecMass2'
        matP = - 1i*k * [matM zeros(size(matM,1)) ; zeros(size(matM,1)) matM];
    case 'PrecDiag'
        matP = matA/diag(diag(matA));
    otherwise
        matP = 1;
end

% Preconditioning (1) — Right
% matA = matA/matP;

% Preconditioning (1) — Left
% matA = matP\matA;
% rhsA = matP\rhsA;
% matP = 1;

% Preconditioning (1) — Symmetric
matP = sqrt(matP);
matA = matP\(matA/matP);
rhsA = matP\rhsA;

% Full system
sysA.matA = matA;
sysA.rhsA = rhsA;

% Compute solution
solA = matA\rhsA;
solA = matP\solA;

% Preconditioning (2)
solA = matP\solA;

% TO IMPROVE IN THE FUTURE
sysA.matP = 1;
sysA.matPinv = 1;

end