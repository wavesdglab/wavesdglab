function [solA, sysA] = computeSolNum1D_DG(mesh, dofm, tau, PREC)

global k BCLeft BCRight

% -------------------------------------------------------------------------
% Quadrature and shape functions
% -------------------------------------------------------------------------

Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);
shapeFunc = functionsShape1D(nodes,dofm.degree);

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

[matElemM, ~, matElemD] = buildMatrixElem1D(dofm.degree);

matM = sparse(dofm.numDof, dofm.numDof);
matD = sparse(dofm.numDof, dofm.numDof);
rhsP = zeros(dofm.numDof,1);
rhsU = zeros(dofm.numDof,1);

for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    [~, ~, ~, ~, souP, souU] = mySol1D(coordGlo);
    
    matLocM = matElemM * length/2;
    matLocD = matElemD;
    
    glo = dofm.locToGlo(e,:);
    matM(glo,glo) = matM(glo,glo) + matLocM;
    matD(glo,glo) = matD(glo,glo) + matLocD;
    rhsP(glo) = rhsP(glo) + (shapeFunc .* souP).' * weights * (length/2);
    rhsU(glo) = rhsU(glo) + (shapeFunc .* souU).' * weights * (length/2);
end

matA = [ -1i*k*matM -matD' ; -matD' -1i*k*matM ];
rhsA = [ rhsP ; rhsU ];

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

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

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

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

sysA.matA = matA;
sysA.rhsA = rhsA;

% Compute solution

solA = sysA.matA\sysA.rhsA;
solA = matP\solA;

end