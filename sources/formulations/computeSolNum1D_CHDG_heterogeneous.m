function [solA, sysA] = computeSolNum1D_CHDG_heterogeneous(mesh, dofm, PREC)

global omega BCLeft BCRight

% -------------------------------------------------------------------------
% Solution at boundaries
% -------------------------------------------------------------------------

% Left
[etaL, ~, cL] = physical_parameters_1D(mesh, 1);
kL = omega / cL;
[solPL, ~, solUL] = mySol1D_heterogeneous(0,kL,etaL);

% Right
[etaR, ~, cR] = physical_parameters_1D(mesh, mesh.numE);
kR = omega / cR;
[solPR, ~, solUR] = mySol1D_heterogeneous(mesh.coordV(mesh.numV),kR,etaR);

[~, ~, c1] = physical_parameters_1D(mesh, 1);
[~, ~, c2] = physical_parameters_1D(mesh, mesh.numE);

coefP = max(c1/c2,c2/c1);    % -> min(c1,c2)/max(c1,c2) < 1

% -------------------------------------------------------------------------
% Quadrature and shape functions
% -------------------------------------------------------------------------

Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeQ = functionsShapeLIN(nodes,dofm.degree);
shapeDerQ = functionsShapeDerLIN(nodes,dofm.degree);
matMelem = shapeQ' * (weights .* shapeQ);
matDelem = shapeQ' * (weights .* shapeDerQ);

% -------------------------------------------------------------------------
% Memory storage
% -------------------------------------------------------------------------

matII = sparse(2*dofm.numDof, 2*dofm.numDof);
matIG = sparse(2*dofm.numDof, 2*mesh.numE);
matGI = sparse(2*mesh.numE, 2*dofm.numDof);
matGG = sparse(1:2*mesh.numE, 1:2*mesh.numE, 1);

matIIinv = sparse(2*dofm.numDof, 2*dofm.numDof);
matPP = sparse(1:2*mesh.numE, 1:2*mesh.numE, 1);

rhsI = zeros(2*dofm.numDof,1);
rhsG = zeros(2*mesh.numE, 1);

for e=1:mesh.numE

    % -------------------------------------------------------------------------
    % Build local systems
    % -------------------------------------------------------------------------

    % Mapping
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;

    % Physical parameters
    [eta, ~, c] = physical_parameters_1D(mesh, e);
    k = omega / c;
    
    % Local RHS vector
    [~, ~, ~, ~, souP, souU] = mySol1D_heterogeneous(coordGlo,k,eta);
    rhsPloc = (shapeQ .* souP).' * weights * (length/2);
    rhsUloc = (shapeQ .* souU).' * weights * (length/2);
    
    % Local matrix
    matMloc = matMelem * length/2;
    matDloc = matDelem;
    matIIloc = [ -1i*k/eta*matMloc -matDloc' ; -matDloc' -1i*k*eta*matMloc ];
    
    % Left local BC
    idIntP = 1;
    idIntU = 1+dofm.numDofPerE;
    if e > 1
        [etaNeigh, ~, ~] = physical_parameters_1D(mesh, e-1);
    else
        etaNeigh = eta;
    end
    matIIloc(idIntP,idIntP) = matIIloc(idIntP,idIntP) + 1/(eta + etaNeigh);          
    matIIloc(idIntP,idIntU) = matIIloc(idIntP,idIntU) - eta/(eta + etaNeigh);          
    matIIloc(idIntU,idIntP) = matIIloc(idIntU,idIntP) - etaNeigh/(eta + etaNeigh);
    matIIloc(idIntU,idIntU) = matIIloc(idIntU,idIntU) + etaNeigh*eta/(eta + etaNeigh);
    
    % Right local BC
    idIntP = 2;
    idIntU = 2+dofm.numDofPerE;
    if e < mesh.numE
        [etaNeigh, ~, ~] = physical_parameters_1D(mesh, e+1);
    else
        etaNeigh = eta;
    end
    matIIloc(idIntP,idIntP) = matIIloc(idIntP,idIntP) + 1/(eta + etaNeigh);
    matIIloc(idIntP,idIntU) = matIIloc(idIntP,idIntU) + eta/(eta + etaNeigh);
    matIIloc(idIntU,idIntP) = matIIloc(idIntU,idIntP) + etaNeigh/(eta + etaNeigh);
    matIIloc(idIntU,idIntU) = matIIloc(idIntU,idIntU) + etaNeigh*eta/(eta + etaNeigh);

    % Assembling
    glo = [dofm.locToGlo(e,:) dofm.locToGlo(e,:)+dofm.numDof]';
    matII(glo,glo) = matIIloc;
    matIIinv(glo,glo) = inv(matIIloc);
    rhsI(glo) = rhsI(glo) + [rhsPloc ; rhsUloc];

%     coefP = max(c,1/c);

    % -------------------------------------------------------------------------
    % Build characteristic variables
    % -------------------------------------------------------------------------

    % Left
    idChar = 2*(e-1)+1;
    idIntP = dofm.locToGlo(e,1);
    idIntU = dofm.locToGlo(e,1) + dofm.numDof;

    if(e > 1)
        [etaNeigh, ~, ~] = physical_parameters_1D(mesh, e-1);
    else
        etaNeigh = eta;
    end

    matIG(idIntP,idChar) = - 1 / (eta + etaNeigh);
    matIG(idIntU,idChar) = - eta / (eta + etaNeigh);

    if(e > 1)
        idIntPneigh = dofm.locToGlo(e-1,2);
        idIntUneigh = dofm.locToGlo(e-1,2) + dofm.numDof;
        matGI(idChar,idIntPneigh) = -1;
        matGI(idChar,idIntUneigh) = -etaNeigh;
    else
        switch BCLeft
            case 'DIR'
                matGI(idChar,idIntP) = +1;
                matGI(idChar,idIntU) = -eta;
                rhsG(idChar) = 2*solPL;
            case 'DIRu'
                matGI(idChar,idIntP) = -1;
                matGI(idChar,idIntU) = eta;
                rhsG(idChar) = 2*solUL;
            case 'ABC'
                matGI(idChar,idIntP) = 0;
                matGI(idChar,idIntU) = 0;
                rhsG(idChar) = solPL + eta * solUL;
            otherwise
                warning('Error - No valid BC has been set on the left.')
        end
    end

    % Right
    idChar = 2*(e-1)+2;
    idIntP = dofm.locToGlo(e,2);
    idIntU = dofm.locToGlo(e,2) + dofm.numDof;

    if(e < mesh.numE)
        [etaNeigh, ~, ~] = physical_parameters_1D(mesh, e+1);
    else
        etaNeigh = eta;
    end

    matIG(idIntP,idChar) = - 1 / (eta + etaNeigh);
    matIG(idIntU,idChar) = eta / (eta + etaNeigh);

    if(e < mesh.numE)
        idIntPneigh = dofm.locToGlo(e+1,1);
        idIntUneigh = dofm.locToGlo(e+1,1) + dofm.numDof;
        matGI(idChar,idIntPneigh) = -1;
        matGI(idChar,idIntUneigh) = etaNeigh;
    else
        switch BCRight
            case 'DIR'
                matGI(idChar,idIntP) = +1;
                matGI(idChar,idIntU) = eta;
                rhsG(idChar) = 2*solPR;
            case 'DIRu'
                matGI(idChar,idIntP) = -1;
                matGI(idChar,idIntU) = - eta;
                rhsG(idChar) = -2*solUR;
            case 'ABC'
                matGI(idChar,idIntP) = 0;
                matGI(idChar,idIntU) = 0;
                rhsG(idChar) = solPR - eta * solUR;
            otherwise
                warning('Error - No valid BC has been set on the left.')
        end
    end

    matPP(2*e-1,2*e-1) = coefP;
    matPP(2*e,2*e) = coefP;
    matPPInv(2*e-1,2*e-1) = inv(coefP);
    matPPInv(2*e,2*e) = inv(coefP);

end

% -------------------------------------------------------------------------
% Build and solve full system
% -------------------------------------------------------------------------

% Build full system
matA = [ matII matIG ; matGI matGG ];
rhsA = [ rhsI ; rhsG ];

% Build reduced system
matS = matGG - matGI*(matIIinv*matIG);
rhsS = rhsG - matGI*(matIIinv*rhsI);

% Build physical system
matPhy = matII - matIG*(matGG\matGI);
rhsPhy = rhsI - matIG*(matGG\rhsG);

% Compute solution
solG = matS\rhsS;
solA = matIIinv*(rhsI-matIG*solG);

% Save system
sysA.matII = matII;
sysA.matIG = matIG;
sysA.matGI = matGI;
sysA.matGG = matGG;
sysA.matIIinv = matIIinv;
sysA.matGGinv = inv(matGG);
sysA.rhsI = rhsI;
sysA.rhsG = rhsG;
sysA.matA = matA;
sysA.rhsA = rhsA;
sysA.matS = matS;
sysA.rhsS = rhsS;
sysA.matPhy = matPhy;
sysA.rhsPhy = rhsPhy;

if PREC == 0
    sysA.matP = 1;
    sysA.matPinv = 1;
else
    sysA.matP = matPP;
    sysA.matPinv = inv(matPP);
end

end
