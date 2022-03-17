function [solFull, matA, rhsA, solRedu, matS, rhsS, matII, matIG, rhsI, matP] ...
    = computeSolNum1D_CG1(mesh, dofm, mySolP, mySolU, mySouP, mySouU)

global k BCLeft BCRight

% Build matrix of the full system

[matM, matK, matD] = buildMatrixGlo1D_CG(mesh, dofm);
matA = [ -1i*k*matM -matD' ; matD -1i*k*matM ];

% Build RHS vector of the full system

rhsP = buildVectorGloRhs1D_CG(mesh, dofm, mySouP);
rhsU = buildVectorGloRhs1D_CG(mesh, dofm, mySouU);
rhsA = [ rhsP ; rhsU ];

% Preconditionning

matP = sparse(1:size(matA,1), 1:size(matA,2), 1);

% Build boundary conditions

if(~(strcmp(BCLeft,'PER') && strcmp(BCRight,'PER')))
    
    switch BCLeft
        case 'DIR'
            rhsA = rhsA - matA(:,1)*mySolP(0);
            rhsA(1) = mySolP(0);
            matA(1,:) = 0;
            matA(:,1) = 0;
            matA(1,1) = 1;
        case 'NEU'
            rhsA(1) = rhsA(1) + mySolU(0);
        case 'ABC'
            rhsA(1) = rhsA(1) + mySolP(0) + mySolU(0);
            matA(1,1) = matA(1,1) + 1;
        otherwise
            warning('Error - No valid BC has been set on the left.')
    end
    
    switch BCRight
        case 'DIR'
            rhsA = rhsA - matA(:,mesh.numV)*mySolP(mesh.coordV(mesh.numV));
            rhsA(mesh.numV) = mySolP(mesh.coordV(mesh.numV));
            matA(mesh.numV,:) = 0;
            matA(:,mesh.numV) = 0;
            matA(mesh.numV,mesh.numV) = 1;
        case 'NEU'
            rhsA(mesh.numV) = rhsA(mesh.numV) - mySolU(mesh.coordV(mesh.numV));
        case 'ABC'
            rhsA(mesh.numV) = rhsA(mesh.numV) + mySolP(mesh.coordV(mesh.numV)) - mySolU(mesh.coordV(mesh.numV));
            matA(mesh.numV,mesh.numV) = matA(mesh.numV,mesh.numV) + 1;
        otherwise
            warning('Error - No valid BC has been set on the right.')
    end
    
end

% Compute solution (full system)

solFull = matA\rhsA;

% Build reduced system

gloGam = [1:dofm.numDofGam (1:dofm.numDofGam)+dofm.numDof];
gloInt = [1:dofm.numDofInt (1:dofm.numDofInt)+dofm.numDof]+dofm.numDofGam;

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

end