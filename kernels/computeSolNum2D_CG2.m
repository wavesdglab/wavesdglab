function [solA, matA, rhsA] = computeSolNum2D_CG2(mesh, dofm)
fprintf('Solver  : Call computeSolNum2D_CG2\n');

global k

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

[matM, matK, ~, ~, rhsA] = buildMatrixGlo2D_CG(mesh, dofm);

matA = matK - k^2*matM;

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

% Quadrature
degreeQ = 4*dofm.degree;
[uQ, weights] = quadratureGaussLIN(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);

% Shape functions
shapeQ = functionsShapeLIN(uQ, dofm.degree);

dofDIR = [];
for edgBnd=1:mesh.numEdgBnd
    edg = mesh.listEdgBnd(edgBnd);
    dof = dofm.locToGloBND(edgBnd,:);
    
    % Normal
    tri = mesh.mapEdgToTri(edg,1);
    fac = mesh.mapEdgToFac(edg,1);
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    normal = getNormalTRI(V1,V2,V3);
    normal = normal(fac,:);
    
    % Mapping
    ver = mesh.mapEdgToVer(edg,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    [xQ, yQ] = locToGloLIN(uQ, V1, V2);
    Jdxdu = norm(V2-V1) * 0.5;  % [ dx/du ]
    
    % Solution function
    [solQ, solDxQ, solDyQ, ~] = mySol(xQ, yQ);
    dirQ = solQ;
    neuQ = normal(1)*solDxQ + normal(2)*solDyQ;
    
    % Orientation
    orientation = ones(dofm.numDofPerLIN,1);
    if(ver(1) > ver(2))
        orientation(3:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    orientation = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation);
    
    % Shape function
    shapeOrQ = shapeQ * orientation;
    
    % Elemental matrix
    matMel = shapeOrQ' * weights * shapeOrQ * Jdxdu;
    rhsDel = shapeOrQ' * weights * dirQ * Jdxdu;
    rhsNel = shapeOrQ' * weights * neuQ * Jdxdu;
    
    % Boundary condition
    switch tagToBC(mesh.tagEdgBnd(edgBnd))
        case 'DIR'
            dofDIR = [dofDIR ; dof];
        case 'NEU'
            rhsA(dof) = rhsA(dof) + rhsNel;
        case 'ABC'
            matA(dof,dof) = matA(dof,dof) - 1i*k * matMel;
            rhsA(dof) = rhsA(dof) + rhsNel - 1i*k * rhsDel;
        otherwise
            warning('Error - No valid BC has been set on the South.')
    end
end

if(~isempty(dofDIR))
    solP = computeSolProjL2_2D_CG(mesh, dofm);
    dofDIR = unique(dofDIR);
    rhsA = rhsA - matA(:,dofDIR)*solP(dofDIR);
    rhsA(dofDIR) = solP(dofDIR);
    matA(dofDIR,:) = 0;
    matA(:,dofDIR) = 0;
    matA(dofDIR,dofDIR) = eye(size(dofDIR,1),size(dofDIR,1));
end

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

solA = matA\rhsA;

end

function BC = tagToBC(tag)
global BCWest BCNorth BCEast BCSouth;
switch tag
    case 1
        BC = BCWest;
    case 2
        BC = BCNorth;
    case 3
        BC = BCEast;
    case 4
        BC = BCSouth;
    otherwise
        warning('Error - No valid BC has been set on the East.')
end
end