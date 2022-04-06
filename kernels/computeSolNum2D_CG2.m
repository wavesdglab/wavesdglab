function [solA, matA, rhsA] = computeSolNum2D_CG2(mesh, dofm)

global k

matA = sparse(dofm.numDofTRI, dofm.numDofTRI);
rhsA = zeros(dofm.numDofTRI, 1);

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

% Quadrature
degreeQ = 4*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);

% Shape functions (f, dfdu, dfdv)
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);
[shapeDuQ, shapeDvQ] = functionsShapeDerTRI(uQ, vQ, dofm.degree);

for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uQ, vQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    Jdudx = inv(Jdxdu);                 % [ du/dx du/dy ; dv/dx dv/dy ]
    detJdxdu = abs(det(Jdxdu));
    
    % RHS function
    [~, ~, ~, rhsQ] = mySol(xQ, yQ);
    
    % Shape functions (dfdx, dfdy)
    shapeDxQ = shapeDuQ * Jdudx(1,1) + shapeDvQ * Jdudx(2,1);
    shapeDyQ = shapeDuQ * Jdudx(1,2) + shapeDvQ * Jdudx(2,2);
    
    % Elemental matrices
    matMel = shapeQ' * weights * shapeQ * detJdxdu;
    matKel = (shapeDxQ' * weights * shapeDxQ + shapeDyQ' * weights * shapeDyQ ) * detJdxdu;
    matAel = matKel - k^2*matMel;
    rhsAel = shapeQ' * weights * rhsQ * detJdxdu;
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matA(dof,dof) = matA(dof,dof) + matAel;
    rhsA(dof) = rhsA(dof) + rhsAel;
    
end

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
    
    % Elemental matrix
    matMel = shapeQ' * weights * shapeQ * Jdxdu;
    rhsDel = shapeQ' * weights * dirQ * Jdxdu;
    rhsNel = shapeQ' * weights * neuQ * Jdxdu;
    
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
    x = mesh.coord(:,1);
    y = mesh.coord(:,2);
    [solP, matP, rhsP] = computeSolProj2D_CG(mesh, dofm);
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