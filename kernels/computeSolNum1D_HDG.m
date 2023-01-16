function [solFull, matA, rhsA, solRedu, matS, rhsS, matIIinv, matIG, rhsI] ...
    = computeSolNum1D_HDG(mesh, dofm, tau)

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

matII = sparse(2*dofm.numDof, 2*dofm.numDof);
matIIinv = sparse(2*dofm.numDof, 2*dofm.numDof);

rhsP = zeros(dofm.numDof,1);
rhsU = zeros(dofm.numDof,1);

for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    
    matMloc = matMelem * length/2;
    matDloc = matDelem;
    
    matIIloc = [ -1i*k*matMloc -matDloc' ; -matDloc' -1i*k*matMloc ];
    
    % Left
    idIntP = 1;
    idIntU = 1+dofm.numDofPerE;
    matIIloc(idIntP,idIntP) = matIIloc(idIntP,idIntP) + tau;
    matIIloc(idIntP,idIntU) = matIIloc(idIntP,idIntU) - 1;
    if((e == 1) && strcmp(BCLeft,'ABC'))
        matIIloc(idIntP,idIntP) = matIIloc(idIntP,idIntP) - tau + 1;
    end
    
    % Right
    idIntP = 2;
    idIntU = 2+dofm.numDofPerE;
    matIIloc(idIntP,idIntP) = matIIloc(idIntP,idIntP) + tau;
    matIIloc(idIntP,idIntU) = matIIloc(idIntP,idIntU) + 1;
    if((e == mesh.numE) && strcmp(BCRight,'ABC'))
        matIIloc(idIntP,idIntP) = matIIloc(idIntP,idIntP) - tau + 1;
    end
    
    idGloP = dofm.locToGlo(e,:);
    idGloU = dofm.locToGlo(e,:) + dofm.numDof;
    matII([idGloP idGloU],[idGloP idGloU]) = matIIloc;
    matIIinv([idGloP idGloU],[idGloP idGloU]) = inv(matIIloc);
    
    % Source term
    [~, ~, ~, ~, souP, souU] = mySol1D(coordGlo);
    
    % Assembling
    glo = dofm.locToGlo(e,:);
    rhsP(glo) = rhsP(glo) + (shapeQ .* souP).' * weights * (length/2);
    rhsU(glo) = rhsU(glo) + (shapeQ .* souU).' * weights * (length/2);
end
rhsI = [ rhsP ; rhsU ];

% Build coupling variables (Pstar)

matGG = sparse(mesh.numV, mesh.numV);
matGI = sparse(mesh.numV, 2*dofm.numDof);
matIG = sparse(2*dofm.numDof, mesh.numV);
rhsG = zeros(mesh.numV, 1);

[solPL, ~, solUL] = mySol1D(0);
[solPR, ~, solUR] = mySol1D(mesh.coordV(mesh.numV));

for e=1:mesh.numE
    
    % Left
    
    idPstar = e;
    idIntP = dofm.locToGlo(e,1);
    idIntU = dofm.locToGlo(e,1) + dofm.numDof;
    matIG(idIntP,idPstar) = matIG(idIntP,idPstar) - tau;
    matIG(idIntU,idPstar) = matIG(idIntU,idPstar) - 1;
    if((e == 1) && strcmp(BCLeft,'ABC'))
        matIG(idIntP,idPstar) = matIG(idIntP,idPstar) + tau - 1;
    end
    
    if(e > 1)
        matGI(idPstar,idIntP)  = matGI(idPstar,idIntP)  + tau;
        matGI(idPstar,idIntU)  = matGI(idPstar,idIntU)  - 1;
        matGG(idPstar,idPstar) = matGG(idPstar,idPstar) - tau;
    else
        switch BCLeft
            case 'DIR'
                matGG(idPstar,idPstar) = matGG(idPstar,idPstar) + 1;
                rhsG(idPstar) = solPL;
            case 'DIRu'
                matGI(idPstar,idIntP)  = matGI(idPstar,idIntP)  + tau;
                matGI(idPstar,idIntU)  = matGI(idPstar,idIntU)  - 1;
                matGG(idPstar,idPstar) = matGG(idPstar,idPstar) - tau;
                rhsG(idPstar) = -solUL;
            case 'ABC'
                matGG(idPstar,idPstar) = matGG(idPstar,idPstar) - 1;
                matGI(idPstar,idIntP)  = matGI(idPstar,idIntP)  + 0.5;
                matGI(idPstar,idIntU)  = matGI(idPstar,idIntU)  - 0.5;
            otherwise
                warning('Error - No valid BC has been set on the left.')
        end
    end
    
    % Right
    
    idPstar = e+1;
    idIntP = dofm.locToGlo(e,2);
    idIntU = dofm.locToGlo(e,2) + dofm.numDof;
    matIG(idIntP,idPstar) = matIG(idIntP,idPstar) - tau;
    matIG(idIntU,idPstar) = matIG(idIntU,idPstar) + 1;
    if((e == mesh.numE) && strcmp(BCRight,'ABC'))
        matIG(idIntP,idPstar) = matIG(idIntP,idPstar) + tau - 1;
    end
    
    if(e < mesh.numE)
        matGI(idPstar,idIntP)  = matGI(idPstar,idIntP)  + tau;
        matGI(idPstar,idIntU)  = matGI(idPstar,idIntU)  + 1;
        matGG(idPstar,idPstar) = matGG(idPstar,idPstar) - tau;
    else
        switch BCRight
            case 'DIR'
                matGG(idPstar,idPstar) = matGG(idPstar,idPstar) + 1;
                rhsG(idPstar) = solPR;
            case 'DIRu'
                matGI(idPstar,idIntP)  = matGI(idPstar,idIntP)  + tau;
                matGI(idPstar,idIntU)  = matGI(idPstar,idIntU)  + 1;
                matGG(idPstar,idPstar) = matGG(idPstar,idPstar) - tau;
                rhsG(idPstar) = solUR;
            case 'ABC'
                matGG(idPstar,idPstar) = matGG(idPstar,idPstar) - 1;
                matGI(idPstar,idIntP)  = matGI(idPstar,idIntP)  + 0.5;
                matGI(idPstar,idIntU)  = matGI(idPstar,idIntU)  + 0.5;
            otherwise
                warning('Error - No valid BC has been set on the right.')
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