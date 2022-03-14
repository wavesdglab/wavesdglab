function [solA, matA, rhsA] = computeSolNumCG2(mesh, dofm)

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

[matM, matK, ~, ~, matS1, matS2, matS3, matS4, dofS1, dofS2, dofS3, dofS4] = buildMatrixGloCG(mesh, dofm);
matA = matK - k^2*matM;
rhsA = matM * solF;

% Build boundary conditions

dofDIR = [];

switch BCWest
    case 'DIR'
        dofDIR = [dofDIR ; dofS1];
    case 'NEU'
        rhsA = rhsA + matS1*(-solDx);
    case 'ABC'
        matA = matA - 1i*k * matS1;
        rhsA = rhsA + matS1*(-solDx - 1i*k*sol);
    otherwise
        warning('Error - No valid BC has been set on the West.')
end

switch BCNorth
    case 'DIR'
        dofDIR = [dofDIR ; dofS2];
    case 'NEU'
        rhsA = rhsA + matS2*(solDy);
    case 'ABC'
        matA = matA - 1i*k * matS2;
        rhsA = rhsA + matS2*(solDy - 1i*k*sol);
    otherwise
        warning('Error - No valid BC has been set on the North.')
end

switch BCEast
    case 'DIR'
        dofDIR = [dofDIR ; dofS3];
    case 'NEU'
        rhsA = rhsA + matS3*(solDx);
    case 'ABC'
        matA = matA - 1i*k * matS3;
        rhsA = rhsA + matS3*(solDx - 1i*k*sol);
    otherwise
        warning('Error - No valid BC has been set on the East.')
end

switch BCSouth
    case 'DIR'
        dofDIR = [dofDIR ; dofS4];
    case 'NEU'
        rhsA = rhsA + matS4*(-solDy);
    case 'ABC'
        matA = matA - 1i*k * matS4;
        rhsA = rhsA + matS4*(-solDy - 1i*k*sol);
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

end