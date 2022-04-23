function [solA, sysA] = computeSolNum2D_UDG1(mesh, dofm, tau)

tic

global k

numDofTRI = dofm.numDofTRI;
numDofFAC = dofm.numDofFAC;
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
matIG    = sparse(3*numDofTRI, numDofFAC);
matGI    = sparse(numDofFAC, 3*numDofTRI);
matGG    = sparse(numDofFAC, numDofFAC);
matIIinv = sparse(3*numDofTRI, 3*numDofTRI);
rhsI     = zeros(3*numDofTRI,1);
rhsG     = zeros(numDofFAC,1);

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
        -1i*k*matMel  -matDXel                            -matDYel                          ;
        -matDXel      -1i*k*matMel                        sparse(numDofPerTRI,numDofPerTRI) ;
        -matDYel       sparse(numDofPerTRI,numDofPerTRI)  -1i*k*matMel                      ];
    
    rhsIel = [
        -1/(1i*k)*vecRHSel    ;
        zeros(numDofPerTRI,1) ;
        zeros(numDofPerTRI,1) ];
    
    % -------------------------------------------------------------------------
    % Surface terms
    % -------------------------------------------------------------------------
    
    matIGel = zeros(3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);
    matGIel = zeros(3*dofm.numDofPerLIN,3*dofm.numDofPerTRI);
    matGGel = zeros(3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);
    rhsGel  = zeros(3*dofm.numDofPerLIN,1);
    
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
        rhsPel = shapeOrQ' * weightsLinQ * solQ * Jdxdu;
        rhsUel = shapeOrQ' * weightsLinQ * solDxQ * Jdxdu / (1i*k);
        rhsVel = shapeOrQ' * weightsLinQ * solDyQ * Jdxdu / (1i*k);
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        % Local ID for interior unknowns and incoming characteristics
        idLocP = 0*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocU = 1*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocV = 2*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocG = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;
        
        % Interior contributions
        matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 0.5*tau           * matMel;
        matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + 0.5     * nx      * matMel;
        matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + 0.5     * ny      * matMel;
        matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 0.5               * matMel;
        
        matIIel(idLocU,idLocP) = matIIel(idLocU,idLocP) + 0.5          * nx * matMel;
        matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + 0.5/tau * nx * nx * matMel;
        matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + 0.5/tau * nx * ny * matMel;
        matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + 0.5/tau      * nx * matMel;
        
        matIIel(idLocV,idLocP) = matIIel(idLocV,idLocP) + 0.5          * ny * matMel;
        matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + 0.5/tau * nx * ny * matMel;
        matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + 0.5/tau * ny * ny * matMel;
        matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + 0.5/tau      * ny * matMel;
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        matGGel(idLocG,idLocG) = matMel;
        
        if (triNeigh > 0)
            
            % Get global ID for exterior unknowns
            idIncG = dofm.locToGloFAC(tri,(1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN);
            dofExt = dofm.locFacNeigh(facNeigh,:);
            idExtP = 0*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = 1*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtV = 2*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            
            matGI(idIncG,idExtP) = matGI(idIncG,idExtP) - tau * matMel;
            matGI(idIncG,idExtU) = matGI(idIncG,idExtU) + nx  * matMel;
            matGI(idIncG,idExtV) = matGI(idIncG,idExtV) + ny  * matMel;
            
        else
            
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            switch tagToBC(mesh.tagEdg(edgGlo))
                case 'DIR'
                    matGIel(idLocG,idLocP) = matGIel(idLocG,idLocP) + tau * matMel;
                    matGIel(idLocG,idLocU) = matGIel(idLocG,idLocU) + nx  * matMel;
                    matGIel(idLocG,idLocV) = matGIel(idLocG,idLocV) + ny  * matMel;
                    rhsGel(idLocG) = rhsGel(idLocG) + 2*tau*rhsPel;
                case 'NEU'
                    matGIel(idLocG,idLocP) = matGIel(idLocG,idLocP) - tau * matMel;
                    matGIel(idLocG,idLocU) = matGIel(idLocG,idLocU) - nx  * matMel;
                    matGIel(idLocG,idLocV) = matGIel(idLocG,idLocV) - ny  * matMel;
                    rhsGel(idLocG) = rhsGel(idLocG) - 2*(nx*rhsUel + ny*rhsVel);
                case 'ABC'
                    
                    matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 0.5*(1-tau)           * matMel;
                    
                    matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + 0.5*(1-1/tau) * nx * nx * matMel;
                    matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + 0.5*(1-1/tau) * nx * ny * matMel;
                    matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + 0.5*(1-1/tau)      * nx * matMel;
                    
                    matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + 0.5*(1-1/tau) * nx * ny * matMel;
                    matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + 0.5*(1-1/tau) * ny * ny * matMel;
                    matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + 0.5*(1-1/tau)      * ny * matMel;
                    
                    rhsGel(idLocG) = rhsGel(idLocG) + (rhsPel - (nx*rhsUel + ny*rhsVel));
                    
                    %rhsGel(idLocG) = rhsGel(idLocG) + (tau*rhsPel - (nx*rhsUel + ny*rhsVel));
                otherwise
                    warning('Error - Bad BC.')
            end
        end
    end
    
    % -------------------------------------------------------------------------
    % Matrix assembling
    % -------------------------------------------------------------------------
    
    dofGloP = 0*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloU = 1*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloV = 2*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloI = [dofGloP dofGloU dofGloV];
    dofGloG = dofm.locToGloFAC(tri,:);
    
    matIIinv(dofGloI,dofGloI) = inv(matIIel);
    matII(dofGloI,dofGloI) = matIIel;
    matIG(dofGloI,dofGloG) = matIGel;
    matGI(dofGloG,dofGloI) = matGIel;
    matGG(dofGloG,dofGloG) = matGGel;
    rhsI(dofGloI)          = rhsIel;
    rhsG(dofGloG)          = rhsGel;
    
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

toc

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