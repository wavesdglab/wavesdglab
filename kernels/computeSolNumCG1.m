function [solA, matA, rhsA] = computeSolNumCG1(mesh, dofm)

global k BCWest BCNorth BCEast BCSouth

% -------------------------------------------------------------------------
% Get data
% -------------------------------------------------------------------------

x = mesh.coord(:,1);
y = mesh.coord(:,2);
[sol, solDx, solDy, solF] = mySol(x,y);

% -------------------------------------------------------------------------
% Build system
% -------------------------------------------------------------------------

[matM, ~, matDX, matDY, matS1, matS2, matS3, matS4, dofS1, dofS2, dofS3, dofS4] = buildMatrixGloCG(mesh, dofm);

matA = [
    -1i*k*matM  -matDX                                 -matDY ;
    matDX'      -1i*k*matM                             sparse(dofm.numDofTRI,dofm.numDofTRI) ;
    matDY'      sparse(dofm.numDofTRI,dofm.numDofTRI)  -1i*k*matM ];
rhsA = -1/(1i*k)*[
    matM*solF ;
    zeros(dofm.numDofTRI,1) ;
    zeros(dofm.numDofTRI,1) ];

% Build boundary conditions

dofTRI = 1:dofm.numDofTRI;
dofDIR = [];

switch BCWest
    case 'DIR'
        dofDIR = [dofDIR ; dofS1];
    case 'NEU'
        rhsA(dofTRI) = rhsA(dofTRI) - matS1 * (-solDx/(1i*k));
    case 'ABC'
        matA(dofTRI,dofTRI) = matA(dofTRI,dofTRI) + matS1;
        rhsA(dofTRI) = rhsA(dofTRI) + matS1 * (sol - (-solDx)/(1i*k));
    otherwise
        warning('Error - No valid BC has been set on the West.')
end

switch BCNorth
    case 'DIR'
        dofDIR = [dofDIR ; dofS2];
    case 'NEU'
        rhsA(dofTRI) = rhsA(dofTRI) - matS2 * (solDy/(1i*k));
    case 'ABC'
        matA(dofTRI,dofTRI) = matA(dofTRI,dofTRI) + matS2;
        rhsA(dofTRI) = rhsA(dofTRI) + matS2 * (sol - (solDy)/(1i*k));
    otherwise
        warning('Error - No valid BC has been set on the North.')
end

switch BCEast
    case 'DIR'
        dofDIR = [dofDIR ; dofS3];
    case 'NEU'
        rhsA(dofTRI) = rhsA(dofTRI) - matS3 * (solDx/(1i*k));
    case 'ABC'
        matA(dofTRI,dofTRI) = matA(dofTRI,dofTRI) + matS3;
        rhsA(dofTRI) = rhsA(dofTRI) + matS3 * (sol - (solDx)/(1i*k));
    otherwise
        warning('Error - No valid BC has been set on the East.')
end

switch BCSouth
    case 'DIR'
        dofDIR = [dofDIR ; dofS4];
    case 'NEU'
        rhsA(dofTRI) = rhsA(dofTRI) - matS4 * (-solDy/(1i*k));
    case 'ABC'
        matA(dofTRI,dofTRI) = matA(dofTRI,dofTRI) + matS4;
        rhsA(dofTRI) = rhsA(dofTRI) + matS4 * (sol - (-solDy)/(1i*k));
    otherwise
        warning('Error - No valid BC has been set on the South.')
end

if(~isempty(dofDIR))
    dofDIR = unique(dofDIR);
    rhsA = rhsA - matA(:,dofDIR)*sol(dofDIR);
    rhsA(dofDIR) = sol(dofDIR);
    matA(dofDIR,:) = 0;
    matA(:,dofDIR) = 0;
    matA(dofDIR,dofDIR) = eye(size(dofDIR,1),size(dofDIR,1));
end

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

solA = matA\rhsA;
solA = solA(dofTRI, 1);

end