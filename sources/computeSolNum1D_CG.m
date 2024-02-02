% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solA, sysA] = computeSolNum1D_CG(mesh, dofm, PREC, shiftPrec)

global k BCLeft BCRight

% -------------------------------------------------------------------------
% Quadrature and shape functions
% -------------------------------------------------------------------------

Q = 16;
[uQ, weightsQ] = quadratureGaussLIN(Q);
shapeQ = functionsShapeLIN(uQ, dofm.degree);
shapeDerQ = functionsShapeDerLIN(uQ, dofm.degree);
matMelem = transpose(shapeQ) * (weightsQ .* shapeQ);
matKelem = transpose(shapeDerQ) * (weightsQ .* shapeDerQ);

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

matM = sparse(dofm.numDof, dofm.numDof);
matK = sparse(dofm.numDof, dofm.numDof);
rhsA = zeros(dofm.numDof,1);

for e=1:mesh.numE
    
    % Mapping
    V1 = mesh.coordV(mesh.listE(e,1));
    V2 = mesh.coordV(mesh.listE(e,2));
    length = abs(V2 - V1);
    xQ = V1*(1-uQ)/2 + V2*(1+uQ)/2;
    
    % Local RHS vector
    [~, ~, ~, rhsVol, ~, ~] = mySol(xQ);
    rhsAloc = (shapeQ .* rhsVol).' * weightsQ * (length/2);
    
    % Local matrices
    matMloc = matMelem * length/2;
    matKloc = matKelem * 2/length;
    
    % Assembling
    glo = dofm.locToGlo(e,:);
    matM(glo,glo) = matM(glo,glo) + matMloc;
    matK(glo,glo) = matK(glo,glo) + matKloc;
    rhsA(glo) = rhsA(glo) + rhsAloc;
end

matA = matK - k^2 * matM;

% -------------------------------------------------------------------------
% Preconditionning
% -------------------------------------------------------------------------

switch PREC
    case 'PrecMass'
        matP = matM;
    case 'PrecDiag'
        matP = diag(diag(matA));
    case 'PrecShiftLap'
        matP = matK - shiftPrec*k^2 * matM;
    otherwise
        matP = sparse(1:size(matA,1), 1:size(matA,2), 1);
end

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

[solL, derL] = mySol(0);
[solR, derR] = mySol(mesh.coordV(mesh.numV));

if(~(strcmp(BCLeft,'PER') && strcmp(BCRight,'PER')))
    
    switch BCLeft
        case 'DIR'
            rhsA = rhsA - matA(:,1)*solL;
            rhsA(1) = solL;
            matA(1,:) = 0;
            matA(:,1) = 0;
            matA(1,1) = 1;
            if(strcmp(PREC,'PrecShiftLap'))
                matP(1,:) = 0;
                matP(:,1) = 0;
                matP(1,1) = 1;
            end
        case 'NEU'
            rhsA(1) = rhsA(1) - derL;
        case 'ABC'
            rhsA(1) = rhsA(1) - derL - 1i*k*solL;
            matA(1,1) = matA(1,1) - 1i*k;
            if(strcmp(PREC,'PrecShiftLap'))
                matP(1,1) = matP(1,1) - 1i*k;
            end
        otherwise
            warning('Error - No valid BC has been set on the left.')
    end
    
    switch BCRight
        case 'DIR'
            rhsA = rhsA - matA(:,mesh.numV)*solR;
            rhsA(mesh.numV) = solR;
            matA(mesh.numV,:) = 0;
            matA(:,mesh.numV) = 0;
            matA(mesh.numV,mesh.numV) = 1;
            if(strcmp(PREC,'PrecShiftLap'))
                matP(mesh.numV,:) = 0;
                matP(:,mesh.numV) = 0;
                matP(mesh.numV,mesh.numV) = 1;
            end
        case 'NEU'
            rhsA(mesh.numV) = rhsA(mesh.numV) + derR;
        case 'ABC'
            rhsA(mesh.numV) = rhsA(mesh.numV) + derR - 1i*k*solR;
            matA(mesh.numV,mesh.numV) = matA(mesh.numV,mesh.numV) - 1i*k;
            if(strcmp(PREC,'PrecShiftLap'))
                matP(mesh.numV,mesh.numV) = matP(mesh.numV,mesh.numV) - 1i*k;
            end
        otherwise
            warning('Error - No valid BC has been set on the right.')
    end
    
end

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

% Preconditioning (1) — Right
% matA = matA/matP;

% Preconditioning (1) — Left
% matA = matP\matA;
% rhsA = matP\rhsA;
% matP = 1;

% Preconditioning (1) — Symmetric
% matP = sqrt(matP);
% matA = matP\(matA/matP);
% rhsA = matP\rhsA;

% Matrix partition
dofG = 1:dofm.numDofGam;
dofI = (dofm.numDofGam+1):(dofm.numDofGam+dofm.numDofInt);
sysA.matII = matA(dofI,dofI);
sysA.matIG = matA(dofI,dofG);
sysA.matGI = matA(dofG,dofI);
sysA.matGG = matA(dofG,dofG);
sysA.rhsI = rhsA(dofI);
sysA.rhsG = rhsA(dofG);

% Full system
sysA.matA = matA;
sysA.rhsA = rhsA;

% Reduced system
sysA.matS = sysA.matGG - sysA.matGI*(sysA.matII\sysA.matIG);
sysA.rhsS = sysA.rhsG - sysA.matGI*(sysA.matII\sysA.rhsI);

% Preconditionning
if (strcmp(PREC,'PrecNone'))
    sysA.matP = 1;
else
    sysA.matP = matP;
end

% Compute solution
solG = sysA.matS\sysA.rhsS;
solI = sysA.matII\(sysA.rhsI-sysA.matIG*solG);
solA = [ solG ; solI ];

% Preconditioning (2)
% solA = matP\solA;

end