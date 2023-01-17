function [solA, sysA] = computeSolNum1D_CG(mesh, dofm, PREC, shiftPrec)

global k BCLeft BCRight

% -------------------------------------------------------------------------
% Quadrature and shape functions
% -------------------------------------------------------------------------

Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeQ = functionsShape1D(nodes,dofm.degree);
shapeDerQ = functionsShapeDer1D(nodes,dofm.degree);
matMelem = shapeQ' * (weights .* shapeQ);
matKelem = shapeDerQ' * (weights .* shapeDerQ);

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

matM = sparse(dofm.numDof, dofm.numDof);
matK = sparse(dofm.numDof, dofm.numDof);
rhsA = zeros(dofm.numDof,1);

for e=1:mesh.numE
    
    % Mapping
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    % Local RHS vector
    [~, ~, ~, rhsVol, ~, ~] = mySol1D(coordGlo);
    rhsAloc = (shapeQ .* rhsVol).' * weights * (length/2);
    
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
        matP = inv(matM);
    case 'PrecDiag'
        matP = matA/diag(diag(matA));
    case 'PrecShiftLap'
        matP = matK - shiftPrec*k^2 * matM;
    otherwise
        matP = sparse(1:size(matA,1), 1:size(matA,2), 1);
end

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

[solL, derL] = mySol1D(0);
[solR, derR] = mySol1D(mesh.coordV(mesh.numV));

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

% Preconditionning
%matA = matA/matP;

%matA = matP\matA;
%rhsA = matP\rhsA;
%matP = 1;

matP = sqrt(matP);
matA = matP\(matA/matP);
rhsA = matP\rhsA;

% Save system
dofG = 1:dofm.numDofGam;
dofI = (1:dofm.numDofInt) + dofm.numDofGam;
sysA.matII = matA(dofI,dofI);
sysA.matIIinv = inv(sysA.matII);
sysA.matIG = matA(dofI,dofG);
sysA.matGI = matA(dofG,dofI);
sysA.matGG = matA(dofG,dofG);
sysA.matA = matA;
sysA.matS = sysA.matGG - sysA.matGI*(sysA.matIIinv*sysA.matIG);
sysA.rhsI = rhsA(dofI);
sysA.rhsG = rhsA(dofG);
sysA.rhsA = rhsA;
sysA.rhsS = sysA.rhsG - sysA.matGI*(sysA.matIIinv*sysA.rhsI);

% Compute solution
solG = sysA.matS\sysA.rhsS;
solI = sysA.matIIinv*(sysA.rhsI-sysA.matIG*solG);
solA = [ solG ; solI ];

solA = matP\solA;

end