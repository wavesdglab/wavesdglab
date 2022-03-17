function [solFull, matA, rhsA, solRedu, matS, rhsS, matIIinv, matIG, rhsI] ...
    = computeSolNum1D_UDG2(mesh, dofm, mySol, mySolDer, mySou, tau)

global k BCLeft BCRight

% Build RHS vector of the system

rhsI = buildVectorGloRhsDG(mesh, dofm, mySou);

% Build matrix of the system

[matElemM, matElemK, matElemD] = buildMatrixElem1D(dofm.degree);
matII = sparse(dofm.numDof, dofm.numDof);
matIIinv = sparse(dofm.numDof, dofm.numDof);
for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    matLocM = matElemM * length/2;
    matLocK = matElemK * 2/length;
    
    matIIloc = matLocK - k^2 * matLocM;
    
    if(e > 1)
        matIIloc(1,1) = matIIloc(1,1) - 1i*tau*k;
    else
        switch BCLeft
            case 'DIR'
                matIIloc(1,:) = 0;
                matIIloc(1,1) = 1;
                rhsI(1) = mySol(0);
            case 'NEU'
                rhsI(1) = rhsI(1) - mySolDer(0);
            case 'ABC'
                matIIloc(1,1) = matIIloc(1,1) - 1i*k;
            otherwise
                warning('Error - No valid BC has been set on the left.')
        end
    end
    
    if(e < mesh.numE)
        matIIloc(2,2) = matIIloc(2,2) - 1i*tau*k;
    else
        switch BCRight
            case 'DIR'
                matIIloc(2,:) = 0;
                matIIloc(2,2) = 1;
                rhsI(dofm.locToGlo(e,2)) = mySol(mesh.coordV(mesh.numV));
            case 'NEU'
                rhsI(dofm.locToGlo(e,2)) = rhsI(dofm.locToGlo(e,2)) + mySolDer(mesh.coordV(mesh.numV));
            case 'ABC'
                matIIloc(2,2) = matIIloc(2,2) - 1i*k;
            otherwise
                warning('Error - No valid BC has been set on the left.')
        end
    end
    
    idGlo = dofm.locToGlo(e,:);
    matII(idGlo,idGlo) = matIIloc;
    matIIinv(idGlo,idGlo) = inv(matIIloc);
end

% Build characteristic variables

matGG = sparse(1:2*mesh.numE, 1:2*mesh.numE, 1);
matGI = sparse(2*mesh.numE, dofm.numDof);
matIG = sparse(dofm.numDof, 2*mesh.numE);
rhsG = zeros(2*mesh.numE, 1);

for e=1:mesh.numE
    
    % Left
    idChar = 2*(e-1)+1;
    idInt = dofm.locToGlo(e,1);
    
    if(e > 1)
        idCharNeigh = idChar-1;
        matIG(idInt,idChar) = -1;
        matGG(idChar,idCharNeigh) = 1;
        matGI(idChar,idInt) = 2i*tau*k;
    end
    
    % Right
    idChar = 2*(e-1)+2;
    idInt = dofm.locToGlo(e,2);
    
    if(e < mesh.numE)
        idCharNeigh = idChar+1;
        matIG(idInt,idChar) = -1;
        matGG(idChar,idCharNeigh) = 1;
        matGI(idChar,idInt) = 2i*tau*k;
    end
    
end

% Build and solve full system

matA = [ matII matIG ; matGI matGG ];
rhsA = [ rhsI ; rhsG ];
solA = matA\rhsA;
solFull = solA(1:dofm.numDof);

% Build and solve reduced system

matS = matGG - matGI*(matIIinv*matIG);
rhsS = rhsG - matGI*(matIIinv*rhsI);
solG = matS\rhsS;
solI = matIIinv*(rhsI-matIG*solG);
solRedu = solI;

end