function [solFull, matA, rhsA, solRedu, matS, rhsS, matIIinv, matIG, rhsI] ...
    = computeSolNum1D_UDG1(mesh, dofm, mySolP, mySolU, mySouP, mySouU, tau)

global k BCLeft BCRight

% Build matrix of the system

[matElemM, matElemK, matElemD] = buildMatrixElem1D(dofm.degree);
matII = sparse(2*dofm.numDof, 2*dofm.numDof);
matIIinv = sparse(2*dofm.numDof, 2*dofm.numDof);
for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    matLocM = matElemM * length/2;
    matLocD = matElemD;
    
    matIIloc = [ -1i*k*matLocM -matLocD' ; -matLocD' -1i*k*matLocM ];
    
    % Left
    idIntP = 1;
    idIntU = 1+dofm.numDofPerE;
    matIIloc(idIntP,idIntP) = matIIloc(idIntP,idIntP) + 0.5*tau;
    matIIloc(idIntP,idIntU) = matIIloc(idIntP,idIntU) - 0.5;
    matIIloc(idIntU,idIntP) = matIIloc(idIntU,idIntP) - 0.5;
    matIIloc(idIntU,idIntU) = matIIloc(idIntU,idIntU) + 0.5/tau;
    if((e == 1) && strcmp(BCLeft,'ABC'))
        matIIloc(idIntP,idIntP) = matIIloc(idIntP,idIntP) - 0.5*tau + 0.5;
        matIIloc(idIntU,idIntU) = matIIloc(idIntU,idIntU) - 0.5/tau + 0.5;
    end
    
    % Right
    idIntP = 2;
    idIntU = 2+dofm.numDofPerE;
    matIIloc(idIntP,idIntP) = matIIloc(idIntP,idIntP) + 0.5*tau;
    matIIloc(idIntP,idIntU) = matIIloc(idIntP,idIntU) + 0.5;
    matIIloc(idIntU,idIntP) = matIIloc(idIntU,idIntP) + 0.5;
    matIIloc(idIntU,idIntU) = matIIloc(idIntU,idIntU) + 0.5/tau;
    if((e == mesh.numE) && strcmp(BCRight,'ABC'))
        matIIloc(idIntP,idIntP) = matIIloc(idIntP,idIntP) - 0.5*tau + 0.5;
        matIIloc(idIntU,idIntU) = matIIloc(idIntU,idIntU) - 0.5/tau + 0.5;
    end
    
    idGloP = dofm.locToGlo(e,:);
    idGloU = dofm.locToGlo(e,:) + dofm.numDof;
    matII([idGloP idGloU],[idGloP idGloU]) = matIIloc;
    matIIinv([idGloP idGloU],[idGloP idGloU]) = inv(matIIloc);
end

% Build RHS vector of the system

rhsP = buildVectorGloRhs1D_DG(mesh, dofm, mySouP);
rhsU = buildVectorGloRhs1D_DG(mesh, dofm, mySouU);
rhsI = [ rhsP ; rhsU ];

% Build characteristic variables

matGG = sparse(1:2*mesh.numE, 1:2*mesh.numE, 1);
matGI = sparse(2*mesh.numE, 2*dofm.numDof);
matIG = sparse(2*dofm.numDof, 2*mesh.numE);
rhsG = zeros(2*mesh.numE, 1);

for e=1:mesh.numE
    
    % Left
    idChar = 2*(e-1)+1;
    idIntP = dofm.locToGlo(e,1);
    idIntU = dofm.locToGlo(e,1) + dofm.numDof;
    matIG(idIntP,idChar) = -0.5;
    matIG(idIntU,idChar) = -0.5/tau;
    
    if(e > 1)
        idIntPneigh = dofm.locToGlo(e-1,2);
        idIntUneigh = dofm.locToGlo(e-1,2) + dofm.numDof;
        matGI(idChar,idIntPneigh) = -tau;
        matGI(idChar,idIntUneigh) = -1;
    else
        switch BCLeft
            case 'DIR'
                matGI(idChar,idIntP) = +tau;
                matGI(idChar,idIntU) = -1;
                rhsG(idChar) = 2*tau*mySolP(0);
            case 'DIRu'
                matGI(idChar,idIntP) = -tau;
                matGI(idChar,idIntU) = +1;
                rhsG(idChar) = 2*mySolU(0);
            case 'ABC'
                % (nothing to do)
            otherwise
                warning('Error - No valid BC has been set on the left.')
        end
    end
    
    % Right
    idChar = 2*(e-1)+2;
    idIntP = dofm.locToGlo(e,2);
    idIntU = dofm.locToGlo(e,2) + dofm.numDof;
    matIG(idIntP,idChar) = -0.5;
    matIG(idIntU,idChar) =  0.5/tau;
    
    if(e < mesh.numE)
        idIntPneigh = dofm.locToGlo(e+1,1);
        idIntUneigh = dofm.locToGlo(e+1,1) + dofm.numDof;
        matGI(idChar,idIntPneigh) = -tau;
        matGI(idChar,idIntUneigh) = +1;
    else
        switch BCRight
            case 'DIR'
                matGI(idChar,idIntP) = +tau;
                matGI(idChar,idIntU) = +1;
                rhsG(idChar) = 2*tau*mySolP(mesh.coordV(mesh.numV));
            case 'DIRu'
                matGI(idChar,idIntP) = -tau;
                matGI(idChar,idIntU) = -1;
                rhsG(idChar) = -2*mySolU(mesh.coordV(mesh.numV));
            case 'ABC'
                % (nothing to do)
            otherwise
                warning('Error - No valid BC has been set on the left.')
        end
    end
    
end

% Build and solve full system

matA = [ matII matIG ; matGI matGG ];
rhsA = [ rhsI ; rhsG ];
solA = matA\rhsA;
solFull = solA(1:2*dofm.numDof);

% Build and solve reduced system


matS = matGG - matGI*(matIIinv*matIG);
rhsS = rhsG - matGI*(matIIinv*rhsI);

%matP = matS/diag(diag(matS));
%matS = matS/matP;

solG = matS\rhsS;

%solG = matP\solG;

solI = matIIinv*(rhsI-matIG*solG);
solRedu = solI;

end