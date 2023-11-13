% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solA, sysA] = computeSolNum2D_CG(mesh, dofm, PREC)

global k

matA = sparse(dofm.numDofTRI,dofm.numDofTRI);
matM = sparse(dofm.numDofTRI,dofm.numDofTRI);
rhsA = zeros(dofm.numDofTRI, 1);


% -------------------------------------------------------------------------
% Quadrature
% -------------------------------------------------------------------------

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);

% Shape functions
shapeLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeDuQ, shapeDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uTriQ, vTriQ, V1, V2, V3); 
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    Jdudx = inv(Jdxdu);                 % [ du/dx du/dy ; dv/dx dv/dy ]
    detJdxdu = abs(det(Jdxdu));

    
    % Orientation
    orientation = ones(dofm.numDofPerTRI,1);
    if(ver(1) > ver(2))
        orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(ver(2) > ver(3))
        orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(ver(3) > ver(1))
        orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    orientation = sparse(1:dofm.numDofPerTRI, 1:dofm.numDofPerTRI, orientation);
    
    % Shape functions (f, dfdx, dfdy) with orientation
    shapeOrQ = shapeTriQ * orientation;
    shapeDxQ = (shapeDuQ * Jdudx(1,1) + shapeDvQ * Jdudx(2,1)) * orientation;
    shapeDyQ = (shapeDuQ * Jdudx(1,2) + shapeDvQ * Jdudx(2,2)) * orientation;
    
    % RHS function
    [~, ~, ~, rhsQ] = mySol(xQ, yQ);
    
    % Elemental matrices
    matMel = shapeOrQ' * (weightsTriQ .* shapeOrQ) * detJdxdu;
    matKel = (shapeDxQ' * (weightsTriQ .* shapeDxQ) + shapeDyQ' * (weightsTriQ .* shapeDyQ) ) * detJdxdu;
    rhsPel = shapeOrQ' * (weightsTriQ .* rhsQ) * detJdxdu;

    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matA(dof,dof) = matA(dof,dof) + matKel - k^2*matMel;
    matM(dof,dof) = matM(dof,dof) + matMel;
    rhsA(dof) = rhsA(dof) + rhsPel;
    % matShiftedLaplacian(dof,dof) = matShiftedLaplacian(dof,dof) + matKel + k^2*matMel;
    
end

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

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
    [xQ, yQ] = locToGloLIN(uLinQ, V1, V2);
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
    
    % Shape function with orientation
    shapeOrQ = shapeLinQ * orientation;
    
    % Elemental matrices/vectors
    matMel = shapeOrQ' * (weightsLinQ .* shapeOrQ) * Jdxdu;
    rhsDel = shapeOrQ' * (weightsLinQ .* dirQ) * Jdxdu;
    rhsNel = shapeOrQ' * (weightsLinQ .* neuQ) * Jdxdu;
    
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
            error('BAD BOUNDARY CONDITION.');
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

% Matrix partition
numDofTRIred = mesh.numVer * dofm.numDofPerVer + mesh.numEdg * dofm.numDofPerEdg;
dofG = 1:numDofTRIred; % noeuds qu'on garde : noeuds des vertex + arètes
dofI = (numDofTRIred+1):dofm.numDofTRI;
sysA.matII = matA(dofI,dofI); % ddl intérieurs 
sysA.matIG = matA(dofI,dofG);
sysA.matGI = matA(dofG,dofI);
sysA.matGG = matA(dofG,dofG);
sysA.rhsI = rhsA(dofI);
sysA.rhsG = rhsA(dofG);
% sysA.matIIinv = inv(sysA.matII);

% Full system
sysA.matA = matA;
sysA.rhsA = rhsA;

% Reduced system
% sysA.matS = sysA.matGG - sysA.matGI*(sysA.matIIinv*sysA.matIG);
% sysA.rhsS = sysA.rhsG - sysA.matGI*(sysA.matIIinv*sysA.rhsI);

% Preconditionning
if (PREC == 1)
    warning('DO NOT USE Pinv : USE P\ INSTEAD');
    % sysA.matP = matM;
    % sysA.matP = matShiftedLaplacian;
    sysA.matP = 1;
    sysA.matPinv = 1;
else
    sysA.matP = 1;
    sysA.matPinv = 1;
end

% Compute solution
solG = sysA.matS\sysA.rhsS;
solI = sysA.matIIinv*(sysA.rhsI-sysA.matIG*solG);
solA = [ solG ; solI ];
% solA = 0;

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
        error('BAD BOUNDARY TAG.')
end
end