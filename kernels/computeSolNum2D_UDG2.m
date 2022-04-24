function [solA, sysA] = computeSolNum2D_UDG2(mesh, dofm, tau)

tic

global k

% Global matrices
matII    = sparse(dofm.numDofTRI, dofm.numDofTRI);
matIG    = sparse(dofm.numDofTRI, dofm.numDofFAC);
matGI    = sparse(dofm.numDofFAC, dofm.numDofTRI);
matGG    = sparse(dofm.numDofFAC, dofm.numDofFAC);
matIIinv = sparse(dofm.numDofTRI, dofm.numDofTRI);
rhsI     = zeros(dofm.numDofTRI,1);
rhsG     = zeros(dofm.numDofFAC,1);

% -------------------------------------------------------------------------
% Quadrature
% -------------------------------------------------------------------------

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);
weightsLinQ = sparse(1:size(weightsLinQ,1), 1:size(weightsLinQ,1), weightsLinQ);
weightsTriQ = sparse(1:size(weightsTriQ,1), 1:size(weightsTriQ,1), weightsTriQ);

% Shape functions
shapeLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeTriDuQ, shapeTriDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

dofDIR = [];

for tri=1:mesh.numTri
    
    % -------------------------------------------------------------------------
    % Volume terms
    % -------------------------------------------------------------------------
    
    % Mapping
    verTri = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(verTri(1),:);
    V2 = mesh.coord(verTri(2),:);
    V3 = mesh.coord(verTri(3),:);
    [xQ, yQ] = locToGloTRI(uTriQ, vTriQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    Jdudx = inv(Jdxdu);                 % [ du/dx du/dy ; dv/dx dv/dy ]
    detJdxdu = abs(det(Jdxdu));
    
    % Orientation
    orientation = ones(dofm.numDofPerTRI,1);
    if(verTri(1) > verTri(2))
        orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(verTri(2) > verTri(3))
        orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(verTri(3) > verTri(1))
        orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    orientation = sparse(1:dofm.numDofPerTRI, 1:dofm.numDofPerTRI, orientation);
    
    % Shape functions (f, dfdx, dfdy) with orientation
    shapeOrQ = shapeTriQ * orientation;
    shapeDxQ = (shapeTriDuQ * Jdudx(1,1) + shapeTriDvQ * Jdudx(2,1)) * orientation;
    shapeDyQ = (shapeTriDuQ * Jdudx(1,2) + shapeTriDvQ * Jdudx(2,2)) * orientation;
    
    % RHS function
    [~, ~, ~, rhsQ] = mySol(xQ, yQ);
    
    % Elemental matrices
    matMel = shapeOrQ' * weightsTriQ * shapeOrQ * detJdxdu;
    matKel = shapeDxQ' * weightsTriQ * shapeDxQ * detJdxdu + shapeDyQ' * weightsTriQ * shapeDyQ * detJdxdu;
    rhsIel = shapeOrQ' * weightsTriQ * rhsQ * detJdxdu;
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matII(dof,dof) = matII(dof,dof) + matKel - k^2*matMel;
    rhsI(dof) = rhsI(dof) + rhsIel;
    
    % -------------------------------------------------------------------------
    % Surface terms
    % -------------------------------------------------------------------------
    
    % Exterior normals
    normal = getNormalTRI(V1,V2,V3);
    
    n1 = mesh.mapTriToVer(tri,:);
    n2 = [n1(2) n1(3) n1(1)]';
    
    % Loop over faces
    for fac = 1:3
        
        % Mapping
        V1 = mesh.coord(n1(fac),:);
        V2 = mesh.coord(n2(fac),:);
        [xQ, yQ] = locToGloLIN(uLinQ, V1, V2);
        Jdxdu = norm(V2-V1) * 0.5;  % [ dx/du ]
        
        % Solution function
        [solQ, solDxQ, solDyQ, ~] = mySol(xQ, yQ);
        dirQ = solQ;
        neuQ = normal(fac,1)*solDxQ + normal(fac,2)*solDyQ;
        
        % Orientation
        orientation = ones(dofm.numDofPerLIN,1);
        if(n1(fac) > n2(fac))
            orientation(3:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        orientation = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation);
        
        % Shape function with orientation
        shapeOrQ = shapeLinQ * orientation;
        
        % Elemental matrices/vectors
        matMel = shapeOrQ' * weightsLinQ * shapeOrQ * Jdxdu;
        rhsDel = shapeOrQ' * weightsLinQ * dirQ * Jdxdu;
        rhsNel = shapeOrQ' * weightsLinQ * neuQ * Jdxdu;
        
        % Local/Global ID for interior unknowns and incoming characteristics
        idLocPint = dofm.locFac(fac,:);
        idLocGint = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;
        idGloPint = dofm.locToGloTRI(tri,idLocPint);
        idGloGinc = dofm.locToGloFAC(tri,idLocGint);
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        matGG(idGloGinc,idGloGinc) = matMel;
        
        if (triNeigh > 0)
            
            matII(idGloPint,idGloPint) = matII(idGloPint,idGloPint) - 1i*tau*k * matMel;
            matIG(idGloPint,idGloGinc) = matIG(idGloPint,idGloGinc) -            matMel;
            
            % Local/Global ID for outgoing characteristics
            idLocGout = (1:dofm.numDofPerLIN) + (facNeigh-1)*dofm.numDofPerLIN;
            idGloGout = dofm.locToGloFAC(triNeigh,idLocGout);
            
            matGG(idGloGinc,idGloGout) = matGG(idGloGinc,idGloGout) +            matMel;
            matGI(idGloGinc,idGloPint) = matGI(idGloGinc,idGloPint) + 2i*tau*k * matMel;
            
        else
            
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            switch tagToBC(mesh.tagEdg(edgGlo))
                case 'DIR'
                    dofDIR = [dofDIR ; idGloPint];
                case 'NEU'
                    rhsI(idGloPint) = rhsI(idGloPint) + rhsNel;
                case 'ABC'
                    matII(idGloPint,idGloPint) = matII(idGloPint,idGloPint) - 1i*k * matMel;
                    rhsI(idGloPint) = rhsI(idGloPint) + rhsNel - 1i*k*rhsDel;
                otherwise
                    warning('Error - Bad BC.')
            end
        end
    end
end

if(~isempty(dofDIR))
    solP = computeSolProjL2_2D_DG(mesh, dofm);
    dofDIR = unique(dofDIR);
    rhsI = rhsI - matII(:,dofDIR)*solP(dofDIR);
    rhsI(dofDIR) = solP(dofDIR);
    matII(dofDIR,:) = 0;
    matII(:,dofDIR) = 0;
    matII(dofDIR,dofDIR) = eye(size(dofDIR,1),size(dofDIR,1));
end

toc

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

tic

matS = matGG - matGI*(matIIinv*matIG);
rhsS = rhsG - matGI*(matIIinv*rhsI);
solG = matS\rhsS;
solI = matIIinv*(rhsI-matIG*solG);
solA = [ solI ; solG ];

toc

sysA.matIIinv = matIIinv;
sysA.matII = matII;
sysA.matIG = matIG;
sysA.matGI = matGI;
sysA.matGG = matGG;
sysA.matA = [ matII matIG ; matGI matGG ];
sysA.matS = matS;
sysA.rhsI = rhsI;
sysA.rhsG = rhsG;
sysA.rhsA = [ rhsI ; rhsG ];
sysA.rhsS = rhsS;

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