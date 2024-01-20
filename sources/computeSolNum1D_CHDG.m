% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solA, sysA] = computeSolNum1D_CHDG(mesh, dofm, tau)

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
% Build local systems
% -------------------------------------------------------------------------

matII = sparse(2*dofm.numDof, 2*dofm.numDof);
matIIinv = sparse(2*dofm.numDof, 2*dofm.numDof);
rhsI = zeros(2*dofm.numDof,1);
for e=1:mesh.numE
    
    % Mapping
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    % Local RHS vector
    [~, ~, ~, ~, souP, souU] = mySol(coordGlo);
    rhsPloc = (shapeQ .* souP).' * weights * (length/2);
    rhsUloc = (shapeQ .* souU).' * weights * (length/2);
    
    % Local matrix
    matMloc = matMelem * length/2;
    matDloc = matDelem;
    matIIloc = [ -1i*k*matMloc -matDloc' ; -matDloc' -1i*k*matMloc ];
    
    % Left local BC
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
    
    % Right local BC
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
    
    % Assembling
    glo = [dofm.locToGlo(e,:) dofm.locToGlo(e,:)+dofm.numDof]';
    matII(glo,glo) = matIIloc;
    matIIinv(glo,glo) = inv(matIIloc);
    rhsI(glo) = rhsI(glo) + [rhsPloc ; rhsUloc];
end

% -------------------------------------------------------------------------
% Build characteristic variables
% -------------------------------------------------------------------------

[solPL, ~, solUL] = mySol(0);
[solPR, ~, solUR] = mySol(mesh.coordV(mesh.numV));

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
                rhsG(idChar) = 2*tau*solPL;
            case 'NEU'
                matGI(idChar,idIntP) = -tau;
                matGI(idChar,idIntU) = +1;
                rhsG(idChar) = 2*solUL;
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
                rhsG(idChar) = 2*tau*solPR;
            case 'NEU'
                matGI(idChar,idIntP) = -tau;
                matGI(idChar,idIntU) = -1;
                rhsG(idChar) = -2*solUR;
            case 'ABC'
                % (nothing to do)
            otherwise
                warning('Error - No valid BC has been set on the left.')
        end
    end
    
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

% TO IMPROVE IN THE FUTURE
sysA.matP = 1;
sysA.matPinv = 1;

end