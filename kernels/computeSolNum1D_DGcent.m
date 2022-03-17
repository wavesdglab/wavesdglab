function [solFull, matA, rhsA, matP] ...
    = computeSolNum1D_DGcent(mesh, dofm, mySolP, mySolU, mySouP, mySouU, PREC)

global k BCLeft BCRight

% Build matrix of the system

[matM, matK, matD] = buildMatrixGlo1D_DG(mesh, dofm);
matA = [ -1i*k*matM -matD' ; -matD' -1i*k*matM ];

% Build RHS vector of the system

rhsP = buildVectorGloRhs1D_DG(mesh, dofm, mySouP);
rhsU = buildVectorGloRhs1D_DG(mesh, dofm, mySouU);
rhsA = [ rhsP ; rhsU ];

% Build interfaces and boundary conditions

for e=1:mesh.numE
    
    % Left
    idIntP = dofm.locToGlo(e,1);
    idIntU = dofm.locToGlo(e,1) + dofm.numDof;
    
    if(e > 1)
        idExtP = dofm.locToGlo(e-1,2);
        idExtU = dofm.locToGlo(e-1,2) + dofm.numDof;
        matA(idIntP,idIntU) = matA(idIntP,idIntU) - 0.5;
        matA(idIntP,idExtU) = matA(idIntP,idExtU) - 0.5;
        matA(idIntU,idIntP) = matA(idIntU,idIntP) - 0.5;
        matA(idIntU,idExtP) = matA(idIntU,idExtP) - 0.5;
    else
        switch BCLeft
            case 'PER'
                idExtP = dofm.locToGlo(mesh.numE,2);
                idExtU = dofm.locToGlo(mesh.numE,2) + dofm.numDof;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) - 0.5;
                matA(idIntP,idExtU) = matA(idIntP,idExtU) - 0.5;
                matA(idIntU,idIntP) = matA(idIntU,idIntP) - 0.5;
                matA(idIntU,idExtP) = matA(idIntU,idExtP) - 0.5;
            case 'DIR'
                matA(idIntP,idIntU) = matA(idIntP,idIntU) - 1;
                rhsA(idIntU) = rhsA(idIntU) + mySolP(0);
            case 'NEU'
                matA(idIntU,idIntP) = matA(idIntU,idIntP) - 1;
                rhsA(idIntP) = rhsA(idIntP) + mySolU(0);
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
        matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5;
        matA(idIntP,idExtU) = matA(idIntP,idExtU) + 0.5;
        matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5;
        matA(idIntU,idExtP) = matA(idIntU,idExtP) + 0.5;
    else
        switch BCRight
            case 'PER'
                idExtP = dofm.locToGlo(1,1);
                idExtU = dofm.locToGlo(1,1) + dofm.numDof;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5;
                matA(idIntP,idExtU) = matA(idIntP,idExtU) + 0.5;
                matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5;
                matA(idIntU,idExtP) = matA(idIntU,idExtP) + 0.5;
            case 'DIR'
                matA(idIntP,idIntU) = matA(idIntP,idIntU) + 1;
                rhsA(idIntU) = rhsA(idIntU) - mySolP(mesh.coordV(mesh.numV));
            case 'NEU'
                matA(idIntU,idIntP) = matA(idIntU,idIntP) + 1;
                rhsA(idIntP) = rhsA(idIntP) - mySolU(mesh.coordV(mesh.numV));
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

% Preconditioning

switch PREC
    case 'PrecMass'
        matP = [matM zeros(size(matM,1)) ; zeros(size(matM,1)) matM];
    case 'PrecDiag'
        matP = matA/diag(diag(matA));
    case 'PrecDissip'
        matP = [ 1i*k*matM -matD' ; -matD' 1i*k*matM ];
    otherwise
        matP = 1;
end

matA = matA/matP;

%matA = matP\matA;
%rhsA = matP\rhsA;
%matP = 1;

%matP = sqrt(matP);
%matA = matP\(matA/matP);
%rhsA = matP\rhsA;

% Compute solution

solFull = matA\rhsA;
solFull = matP\solFull;

end