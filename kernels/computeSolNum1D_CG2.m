function [solFull, matA, rhsA, solRedu, matS, rhsS, matII, matIG, rhsI, matP] ...
    = computeSolNum1D_CG2(mesh, dofm, mySol, mySolDer, mySou, PREC, alphaPrec)

global k BCLeft BCRight

% Build matrix of the full system

[matM, matK, ~] = buildMatrixGlo1D_CG(mesh, dofm);
matA = matK - k^2 * matM;

% Build RHS vector of the full system

rhsA = buildVectorGloRhs1D_CG(mesh, dofm, mySou);

% Preconditionning

switch PREC
    case 'PrecMass'
        matP = inv(matM);
    case 'PrecDiag'
        matP = matA/diag(diag(matA));
    case 'PrecShiftLap'
        matP = matK - alphaPrec*k^2 * matM;
    otherwise
        matP = sparse(1:size(matA,1), 1:size(matA,2), 1);
end

% Build boundary conditions

if(~(strcmp(BCLeft,'PER') && strcmp(BCRight,'PER')))
    
    switch BCLeft
        case 'DIR'
            rhsA = rhsA - matA(:,1)*mySol(0);
            rhsA(1) = mySol(0);
            matA(1,:) = 0;
            matA(:,1) = 0;
            matA(1,1) = 1;
            if(strcmp(PREC,'PrecShiftLap'))
                matP(1,:) = 0;
                matP(:,1) = 0;
                matP(1,1) = 1;
            end
        case 'NEU'
            rhsA(1) = rhsA(1) - mySolDer(0);
        case 'ABC'
            rhsA(1) = rhsA(1) - mySolDer(0) - 1i*k*mySol(0);
            matA(1,1) = matA(1,1) - 1i*k;
            if(strcmp(PREC,'PrecShiftLap'))
                matP(1,1) = matP(1,1) - 1i*k;
            end
        otherwise
            warning('Error - No valid BC has been set on the left.')
    end
    
    switch BCRight
        case 'DIR'
            rhsA = rhsA - matA(:,mesh.numV)*mySol(mesh.coordV(mesh.numV));
            rhsA(mesh.numV) = mySol(mesh.coordV(mesh.numV));
            matA(mesh.numV,:) = 0;
            matA(:,mesh.numV) = 0;
            matA(mesh.numV,mesh.numV) = 1;
            if(strcmp(PREC,'PrecShiftLap'))
                matP(mesh.numV,:) = 0;
                matP(:,mesh.numV) = 0;
                matP(mesh.numV,mesh.numV) = 1;
            end
        case 'NEU'
            rhsA(mesh.numV) = rhsA(mesh.numV) + mySolDer(mesh.coordV(mesh.numV));
        case 'ABC'
            rhsA(mesh.numV) = rhsA(mesh.numV) + mySolDer(mesh.coordV(mesh.numV)) - 1i*k*mySol(mesh.coordV(mesh.numV));
            matA(mesh.numV,mesh.numV) = matA(mesh.numV,mesh.numV) - 1i*k;
            if(strcmp(PREC,'PrecShiftLap'))
                matP(mesh.numV,mesh.numV) = matP(mesh.numV,mesh.numV) - 1i*k;
            end
        otherwise
            warning('Error - No valid BC has been set on the right.')
    end
    
end

% Compute solution (full system)

matA = matA/matP;

%matA = matP\matA;
%rhsA = matP\rhsA;
%matP = 1;

%matP = sqrt(matP);
%matA = matP\(matA/matP);
%rhsA = matP\rhsA;

solFull = matA\rhsA;
solFull = matP\solFull;

% Build reduced system

gloGam = 1:dofm.numDofGam;
gloInt = (1:dofm.numDofInt) + dofm.numDofGam;

matGG = matA(gloGam,gloGam);
matGI = matA(gloGam,gloInt);
matIG = matA(gloInt,gloGam);
matII = matA(gloInt,gloInt);
rhsG = rhsA(gloGam);
rhsI = rhsA(gloInt);

matS = matGG - matGI*(matII\matIG);
rhsS = rhsG - matGI*(matII\rhsI);

% Compute solution (reduced system)

solG = matS\rhsS;
solI = matII\(rhsI-matIG*solG);
solRedu = [ solG ; solI ];
solRedu = matP\solRedu;

end