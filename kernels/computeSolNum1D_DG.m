function [solFull, matA, rhsA, matP] ...
    = computeSolNum1D_DG(mesh, dofm, mySolP, mySolU, mySouP, mySouU, tau, PREC)

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
        matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau;
        matA(idIntP,idIntU) = matA(idIntP,idIntU) - 0.5;
        matA(idIntP,idExtP) = matA(idIntP,idExtP) - 0.5*tau;
        matA(idIntP,idExtU) = matA(idIntP,idExtU) - 0.5;
        matA(idIntU,idIntP) = matA(idIntU,idIntP) - 0.5;
        matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau;
        matA(idIntU,idExtP) = matA(idIntU,idExtP) - 0.5;
        matA(idIntU,idExtU) = matA(idIntU,idExtU) - 0.5/tau;
    else
        switch BCLeft
            case 'PER'
                idExtP = dofm.locToGlo(mesh.numE,2);
                idExtU = dofm.locToGlo(mesh.numE,2) + dofm.numDof;
                matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) - 0.5;
                matA(idIntP,idExtP) = matA(idIntP,idExtP) - 0.5*tau;
                matA(idIntP,idExtU) = matA(idIntP,idExtU) - 0.5;
                matA(idIntU,idIntP) = matA(idIntU,idIntP) - 0.5;
                matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau;
                matA(idIntU,idExtP) = matA(idIntU,idExtP) - 0.5;
                matA(idIntU,idExtU) = matA(idIntU,idExtU) - 0.5/tau;
            case 'DIR'
                matA(idIntP,idIntP) = matA(idIntP,idIntP) + tau;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) - 1;
                rhsA(idIntP) = rhsA(idIntP) + mySolP(0)*tau;
                rhsA(idIntU) = rhsA(idIntU) + mySolP(0);
            case 'DIRu'
                matA(idIntU,idIntP) = matA(idIntU,idIntP) - 1;
                matA(idIntU,idIntU) = matA(idIntU,idIntU) + 1/tau;
                rhsA(idIntP) = rhsA(idIntP) + mySolU(0);
                rhsA(idIntU) = rhsA(idIntU) + mySolU(0)/tau;
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
        matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau;
        matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5;
        matA(idIntP,idExtP) = matA(idIntP,idExtP) - 0.5*tau;
        matA(idIntP,idExtU) = matA(idIntP,idExtU) + 0.5;
        matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5;
        matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau;
        matA(idIntU,idExtP) = matA(idIntU,idExtP) + 0.5;
        matA(idIntU,idExtU) = matA(idIntU,idExtU) - 0.5/tau;
    else
        switch BCRight
            case 'PER'
                idExtP = dofm.locToGlo(1,1);
                idExtU = dofm.locToGlo(1,1) + dofm.numDof;
                matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5;
                matA(idIntP,idExtP) = matA(idIntP,idExtP) - 0.5*tau;
                matA(idIntP,idExtU) = matA(idIntP,idExtU) + 0.5;
                matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5;
                matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau;
                matA(idIntU,idExtP) = matA(idIntU,idExtP) + 0.5;
                matA(idIntU,idExtU) = matA(idIntU,idExtU) - 0.5/tau;
            case 'DIR'
                matA(idIntP,idIntP) = matA(idIntP,idIntP) + tau;
                matA(idIntP,idIntU) = matA(idIntP,idIntU) + 1;
                rhsA(idIntP) = rhsA(idIntP) + mySolP(mesh.coordV(mesh.numV))*tau;
                rhsA(idIntU) = rhsA(idIntU) - mySolP(mesh.coordV(mesh.numV));
            case 'DIRu'
                matA(idIntU,idIntP) = matA(idIntU,idIntP) + 1;
                matA(idIntU,idIntU) = matA(idIntU,idIntU) + 1/tau;
                rhsA(idIntP) = rhsA(idIntP) - mySolU(mesh.coordV(mesh.numV));
                rhsA(idIntU) = rhsA(idIntU) + mySolU(mesh.coordV(mesh.numV))/tau;
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
    case 'PrecMass2'
        matP = - 1i*k * [matM zeros(size(matM,1)) ; zeros(size(matM,1)) matM];
    case 'PrecDiag'
        matP = matA/diag(diag(matA));
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