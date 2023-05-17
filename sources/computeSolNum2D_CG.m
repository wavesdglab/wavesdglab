% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solA, sysA] = computeSolNum2D_CG(mesh, dofm, PREC)

global k
global Rexiy
global L_PML
global L

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

    % if(abs(V1(1))> 1 || abs(V1(2))> 1 || abs(V2(1))> 1 || abs(V2(2))> 1 || abs(V3(1))> 1 || abs(V3(2))> 1)
    %     disp('PML')
    % end

    % if(abs(xQ(1))> 1 || abs(xQ(2))> 1 || abs(xQ(3))> 1 || abs(yQ(1))> 1 || abs(yQ(2))> 1 || abs(yQ(3))> 1)
    %     disp('PML')
    % end

    % sigma_x = zeros(size(xQ)) + ((L <= abs(xQ)) .* (abs(xQ) <= L + L_PML) ./ (L + L_PML-abs(xQ))) ;
    % sigma_y = zeros(size(yQ)) + ((L <= abs(yQ)) .* (abs(yQ) <= L + L_PML) ./ (L + L_PML-abs(yQ))) ;
    % sigma_x = 0;
    % sigma_y = 0;
    % gamma_x = ones(size(xQ)) + 1i*sigma_x/k .* (L <= abs(xQ)) .* (abs(xQ) <= L + L_PML);
    % gamma_y = ones(size(yQ)) + 1i*sigma_y/k .* (L <= abs(yQ)) .* (abs(yQ) <= L + L_PML);
    gamma_x = 1;
    gamma_y = 1;
    alpha_PML = gamma_x .* gamma_y;

    
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
    matMel = shapeOrQ' * (weightsTriQ .* shapeOrQ .* alpha_PML) * detJdxdu;
    matKel = (shapeDxQ' * (weightsTriQ .* shapeDxQ .* gamma_y ./ gamma_x) + shapeDyQ' * (weightsTriQ .* shapeDyQ .* gamma_x ./ gamma_y) ) * detJdxdu;
    rhsPel = shapeOrQ' * (weightsTriQ .* rhsQ) * detJdxdu;
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matA(dof,dof) = matA(dof,dof) + matKel - k^2*matMel;
    matM(dof,dof) = matM(dof,dof) + matMel;
    rhsA(dof) = rhsA(dof) + rhsPel;
    
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
            % If DIR homogene, comment this line :
            dofDIR = [dofDIR ; dof];
            % For DIR homogene :
            % matA(dof,:) = 0;
            % matA(dof,dof) = eye(length(dof),length(dof));
            % rhsA(dof) = 0;
        case 'NEU'
            rhsA(dof) = rhsA(dof) + rhsNel;
        case 'ABC'
            matA(dof,dof) = matA(dof,dof) - 1i*k * matMel;
            rhsA(dof) = rhsA(dof) + rhsNel - 1i*k * rhsDel;
        otherwise
            error('BAD BOUNDARY CONDITION.');
    end
end

% If DIR homogene, comment these lines :
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
dofG = 1:numDofTRIred;
dofI = (numDofTRIred+1):dofm.numDofTRI;
sysA.matII = matA(dofI,dofI);
sysA.matIG = matA(dofI,dofG);
sysA.matGI = matA(dofG,dofI);
sysA.matGG = matA(dofG,dofG);
sysA.rhsI = rhsA(dofI);
sysA.rhsG = rhsA(dofG);
sysA.matIIinv = inv(sysA.matII);

% Full system
sysA.matA = matA;
sysA.rhsA = rhsA;

% Reduced system
sysA.matS = sysA.matGG - sysA.matGI*(sysA.matIIinv*sysA.matIG);
sysA.rhsS = sysA.rhsG - sysA.matGI*(sysA.matIIinv*sysA.rhsI);

% Preconditionning
if (PREC == 1)
    % warning('NO PRECONDITIONNING TECHNIQUE CODED YET FOR CG.')
    sysA.matP = matM;
    sysA.matPinv = inv(matM);
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
global BCWest BCNorth BCEast BCSouth BCCircle;
switch tag
    case 1
        BC = BCWest;
    case 2
        BC = BCNorth;
        % BC = BCCircle;
    case 3
        BC = BCEast;
    case 4
        BC = BCSouth;
    case 5
        BC = BCCircle;
    otherwise
        error('BAD BOUNDARY TAG.')
end
end