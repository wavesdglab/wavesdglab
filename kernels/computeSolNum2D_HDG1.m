function [solA, sysA] = computeSolNum2D_HDG1(mesh, dofm, tau)

tic

global k

numDofTRI = dofm.numDofTRI;
numDofLIN = dofm.numDofLIN;
numDofPerTRI = dofm.numDofPerTRI;

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

% Global matrices
matII    = sparse(3*numDofTRI, 3*numDofTRI);
matIG    = sparse(3*numDofTRI, numDofLIN);
matGI    = sparse(numDofLIN, 3*numDofTRI);
matGG    = sparse(numDofLIN, numDofLIN);
matIIinv = sparse(3*numDofTRI, 3*numDofTRI);
rhsI     = zeros(3*numDofTRI,1);
rhsG     = zeros(numDofLIN,1);

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
    matDXel = shapeDxQ' * weightsTriQ * shapeOrQ * detJdxdu;
    matDYel = shapeDyQ' * weightsTriQ * shapeOrQ * detJdxdu;
    vecRHSel = shapeOrQ' * weightsTriQ * rhsQ * detJdxdu;
    
    matIIel = [
        -1i*k*matMel  -matDXel                           -matDYel                          ;
        -matDXel      -1i*k*matMel                       sparse(numDofPerTRI,numDofPerTRI) ;
        -matDYel      sparse(numDofPerTRI,numDofPerTRI)  -1i*k*matMel                      ];
    
    rhsIel = [
        -1/(1i*k)*vecRHSel    ;
        zeros(numDofPerTRI,1) ;
        zeros(numDofPerTRI,1) ];
    
    % Global ID
    dofGloP = 0*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloU = 1*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloV = 2*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloI = [dofGloP dofGloU dofGloV];
    
    % -------------------------------------------------------------------------
    % Surface terms
    % -------------------------------------------------------------------------
    
    % Exterior normals
    verTri = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(verTri(1),:);
    V2 = mesh.coord(verTri(2),:);
    V3 = mesh.coord(verTri(3),:);
    normal = getNormalTRI(V1,V2,V3);
    
    n1 = mesh.mapTriToVer(tri,:);
    n2 = [n1(2) n1(3) n1(1)]';
    
    % Loop over faces
    for fac = 1:3
        
        matIGel = zeros(3*dofm.numDofPerTRI,dofm.numDofPerLIN);
        matGIel = zeros(dofm.numDofPerLIN,3*dofm.numDofPerTRI);
        matGGel = zeros(dofm.numDofPerLIN,dofm.numDofPerLIN);
        rhsGel = zeros(dofm.numDofPerLIN,1);
        
        % Mapping
        V1 = mesh.coord(n1(fac),:);
        V2 = mesh.coord(n2(fac),:);
        [xQ, yQ] = locToGloLIN(uLinQ, V1, V2);
        Jdxdu = norm(V2-V1) * 0.5;  % [ dx/du ]
        
        % Solution function
        [solQ, solDxQ, solDyQ, ~] = mySol(xQ, yQ);
        
        % Orientation
        orientation = ones(dofm.numDofPerLIN,1);
        if(n1(fac) > n2(fac))
            orientation(3:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        orientation = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation);
        
        % Shape function
        shapeOrQ = shapeLinQ * orientation;
        
        % Elemental matrices
        matMel = shapeOrQ' * weightsLinQ * shapeOrQ * Jdxdu;
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        % Local ID for interior unknowns and incoming characteristics
        idLocP = 0*numDofPerTRI + dofm.locFac(fac,:);
        idLocU = 1*numDofPerTRI + dofm.locFac(fac,:);
        idLocV = 2*numDofPerTRI + dofm.locFac(fac,:);
        
        % Surface terms for the volume fields
        matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + tau * matMel;
        matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + nx  * matMel;
        matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + ny  * matMel;
        matIGel(idLocP,:) = -tau * matMel;
        matIGel(idLocU,:) =  nx  * matMel;
        matIGel(idLocV,:) =  ny  * matMel;
        
        % Surface terms for the surface field
        triNeigh = mesh.mapTriToTri(tri,fac);
        if (triNeigh > 0)
            
            matGIel(:,idLocP) = tau * matMel;
            matGIel(:,idLocU) = nx  * matMel;
            matGIel(:,idLocV) = ny  * matMel;
            matGGel = -tau * matMel;
            
        else
            
            rhsPel = shapeOrQ' * weightsLinQ * solQ * Jdxdu;
            rhsUel = shapeOrQ' * weightsLinQ * solDxQ * Jdxdu / (1i*k);
            rhsVel = shapeOrQ' * weightsLinQ * solDyQ * Jdxdu / (1i*k);
            
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            switch tagToBC(mesh.tagEdg(edgGlo))
                case 'DIR'
                    matGGel = matMel;
                    rhsGel = rhsPel;
                case 'NEU'
                    matGIel(:,idLocP) = tau * matMel;
                    matGIel(:,idLocU) = nx  * matMel;
                    matGIel(:,idLocV) = ny  * matMel;
                    matGGel = -tau * matMel;
                    rhsGel = nx*rhsUel + ny*rhsVel;
                case 'ABC'
                    matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + (1-tau) * matMel;
                    matIGel(idLocP,:) = matIGel(idLocP,:) - (1-tau) * matMel;
                    %
                    matGIel(:,idLocP) =       matMel;
                    matGIel(:,idLocU) = nx  * matMel;
                    matGIel(:,idLocV) = ny  * matMel;
                    matGGel = -2 * matMel;
                    rhsGel  = nx*rhsUel + ny*rhsVel - rhsPel;
                otherwise
                    warning('Error - Bad BC.');
            end
        end
        
        % Global ID for edge unknowns
        edgGlo = abs(mesh.mapTriToEdg(tri,fac));
        dofGloG = dofm.locToGloLIN(edgGlo,:);
        if(mesh.mapTriToEdg(tri,fac) < 0)
            tmp = dofGloG;
            dofGloG(1) = tmp(2);
            dofGloG(2) = tmp(1);
        end
        
        % -------------------------------------------------------------------------
        % Matrix assembling
        % -------------------------------------------------------------------------
        
        matIG(dofGloI,dofGloG) = matIG(dofGloI,dofGloG) + matIGel;
        matGI(dofGloG,dofGloI) = matGI(dofGloG,dofGloI) + matGIel;
        matGG(dofGloG,dofGloG) = matGG(dofGloG,dofGloG) + matGGel;
        rhsG(dofGloG) = rhsG(dofGloG) + rhsGel;
        
    end
    
    % -------------------------------------------------------------------------
    % Matrix assembling
    % -------------------------------------------------------------------------
    
    matIIinv(dofGloI,dofGloI) = inv(matIIel);
    matII(dofGloI,dofGloI)    = matIIel;
    rhsI(dofGloI)             = rhsIel;
    
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
        warning('Error - No valid BC has been set.')
end
end