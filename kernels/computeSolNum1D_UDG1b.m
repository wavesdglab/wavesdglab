function [solFull, matA, rhsA, solRedu, matS, rhsS, matIIinv, matIG, rhsI] ...
    = computeSolNum1D_UDG1b(mesh, dofm, tau)

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
    
    idGloP = dofm.locToGlo(e,:);
    idGloU = dofm.locToGlo(e,:) + dofm.numDof;
    matII([idGloP idGloU],[idGloP idGloU]) = matIIloc;
    matIIinv([idGloP idGloU],[idGloP idGloU]) = inv(matIIloc);
end

% Build RHS vector of the system

Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeFunc = functionsShape1D(nodes,dofm.degree);

rhsP = zeros(dofm.numDof,1);
rhsU = zeros(dofm.numDof,1);
for e=1:mesh.numE
    
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    [~, ~, ~, ~, souP, souU] = mySol1D(coordGlo);
    
    glo = dofm.locToGlo(e,:);
    rhsP(glo) = rhsP(glo) + (shapeFunc .* souP).' * weights * (length/2);
    rhsU(glo) = rhsU(glo) + (shapeFunc .* souU).' * weights * (length/2);
end
rhsA = [ rhsP ; rhsU ];

% Build characteristic variables

matGG = sparse(1:2*mesh.numV, 1:2*mesh.numV, 1);
matGI = sparse(2*mesh.numV, 2*dofm.numDof);
matIG = sparse(2*dofm.numDof, 2*mesh.numV);
rhsG = zeros(2*mesh.numV, 1);

for e=1:mesh.numE
    
    % Left
    idIntP = dofm.locToGlo(e,1);
    idIntU = dofm.locToGlo(e,1) + dofm.numDof;
    idCharInt = 2*(e-1)+2;
    idCharExt = 2*(e-1)+1;
    matIG(idIntP,idCharInt) =  0.5;
    matIG(idIntP,idCharExt) = -0.5;
    if((e == 1) && strcmp(BCLeft,'ABC'))
        matIG(idIntU,idCharInt) = -0.5;
        matIG(idIntU,idCharExt) = -0.5;
    else
        matIG(idIntU,idCharInt) = -0.5/tau;
        matIG(idIntU,idCharExt) = -0.5/tau;
    end
    
    % Left
    idIntP = dofm.locToGlo(e,2);
    idIntU = dofm.locToGlo(e,2) + dofm.numDof;
    idCharInt = 2*(e-1)+3;
    idCharExt = 2*(e-1)+4;
    matIG(idIntP,idCharInt) =  0.5;
    matIG(idIntP,idCharExt) = -0.5;
    if((e == mesh.numE) && strcmp(BCRight,'ABC'))
        matIG(idIntU,idCharInt) =  0.5;
        matIG(idIntU,idCharExt) =  0.5;
    else
        matIG(idIntU,idCharInt) =  0.5/tau;
        matIG(idIntU,idCharExt) =  0.5/tau;
    end
    
end

for v=1:mesh.numV
    
    % Left
    idChar = 2*(v-1)+1;
    if(v > 1)
        idIntPneigh = dofm.locToGlo(v-1,2);
        idIntUneigh = dofm.locToGlo(v-1,2) + dofm.numDof;
        if((v == mesh.numV) && strcmp(BCRight,'ABC'))
            matGI(idChar,idIntPneigh) = -1;
            matGI(idChar,idIntUneigh) = -1;
        else
            matGI(idChar,idIntPneigh) = -tau;
            matGI(idChar,idIntUneigh) = -1;
        end
    else
        idIntP = dofm.locToGlo(v,1);
        idIntU = dofm.locToGlo(v,1) + dofm.numDof;
        switch BCLeft
            case 'DIR'
                matGI(idChar,idIntP) = +tau;
                matGI(idChar,idIntU) = -1;
                rhsG(idChar) = 2*tau*mySolP(0);
            case 'NEU'
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
    idChar = 2*(v-1)+2;
    if(v < mesh.numV)
        idIntPneigh = dofm.locToGlo(v,1);
        idIntUneigh = dofm.locToGlo(v,1) + dofm.numDof;
        if((v == 1) && strcmp(BCLeft,'ABC'))
            matGI(idChar,idIntPneigh) = -1;
            matGI(idChar,idIntUneigh) = +1;
        else
            matGI(idChar,idIntPneigh) = -tau;
            matGI(idChar,idIntUneigh) = +1;
        end
    else
        idIntP = dofm.locToGlo(v-1,2);
        idIntU = dofm.locToGlo(v-1,2) + dofm.numDof;
        switch BCRight
            case 'DIR'
                matGI(idChar,idIntP) = +tau;
                matGI(idChar,idIntU) = +1;
                rhsG(idChar) = 2*tau*mySolP(mesh.coordV(mesh.numV));
            case 'NEU'
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
solG = matS\rhsS;
solI = matIIinv*(rhsI-matIG*solG);
solRedu = solI;

end