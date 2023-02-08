% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solA, sysA] = computeSolNum1D_CHDGb(mesh, dofm, tau)

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
    [~, ~, ~, ~, souP, souU] = mySol1D(coordGlo);
    rhsPloc = (shapeQ .* souP).' * weights * (length/2);
    rhsUloc = (shapeQ .* souU).' * weights * (length/2);
    
    % Local matrix
    matMloc = matMelem * length/2;
    matDloc = matDelem;
    matIIloc = [ -1i*k*matMloc -matDloc' ; -matDloc' -1i*k*matMloc ];
    
    % Assembling
    glo = [dofm.locToGlo(e,:) dofm.locToGlo(e,:)+dofm.numDof];
    matII(glo,glo) = matIIloc;
    matIIinv(glo,glo) = inv(matIIloc);
    rhsI(glo) = rhsI(glo) + [rhsPloc ; rhsUloc];
end

% -------------------------------------------------------------------------
% Build characteristic variables
% -------------------------------------------------------------------------

[solPL, ~, solUL] = mySol1D(0);
[solPR, ~, solUR] = mySol1D(mesh.coordV(mesh.numV));

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