% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solI, sysA, condLoc] = computeSolNum2D_CHDG(mesh, dofm, tau, BASIS, PREC)

global k edgTagToBC
global LdomX LdomY LpmlX LpmlY

numDofTRI = dofm.numDofTRI;
numDofFAC = dofm.numDofFAC;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);

% Shape functions and derivatives (reference space)
shapePhyLinQ = functionsShapeLIN(uLinQ, dofm.degree);        % For physical variables (LIN)
shapePhyTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree); % For physical variables (TRI)
shapeAuxLinQ = functionsLegendre(uLinQ, dofm.degree);        % For auxiliary variables (LIN)
[shapeTriDuQ, shapeTriDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% Global matrices
matIIx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIGx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIGy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIGv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matGIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGGx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matIIvInv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matGGvInv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);

% Global RHS vectors
rhsI = zeros(3*numDofTRI,1);
rhsG = zeros(numDofFAC,1);

condLoc = zeros(mesh.numTri,1);

for tri=1:mesh.numTri
    
    % ---------------------------------------------------------------------
    % Volume terms
    % ---------------------------------------------------------------------
    
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
    
    % Shape functions and derivatives with orientation (physical space)
    shapePhyQ = shapePhyTriQ * orientation;
    shapeDxQ = (shapeTriDuQ * Jdudx(1,1) + shapeTriDvQ * Jdudx(2,1)) * orientation;
    shapeDyQ = (shapeTriDuQ * Jdudx(1,2) + shapeTriDvQ * Jdudx(2,2)) * orientation;
    
    % RHS function
    [~, ~, ~, rhsQ] = mySol(xQ, yQ);
    
    % Elemental matrices/vectors
    weightsQ = weightsTriQ .* detJdxdu;
    matMel = transpose(shapePhyQ) * (weightsQ .* shapePhyQ);
    matDXel = transpose(shapeDxQ) * (weightsQ .* shapePhyQ);
    matDYel = transpose(shapeDyQ) * (weightsQ .* shapePhyQ);
    rhsPel = transpose(shapePhyQ) * (weightsQ .* rhsQ);
    
    matIIel = [
        -1i*k*matMel  -matDXel                          -matDYel                         ;
        -matDXel      -1i*k*matMel                      zeros(numDofPerTRI,numDofPerTRI) ;
        -matDYel      zeros(numDofPerTRI,numDofPerTRI)  -1i*k*matMel                     ];
    
    % PML stretching
    if(~isempty(LdomX) && ~isempty(LdomY))
        sigmaPmlX = (LdomX <= abs(xQ))./(LdomX+LpmlX-abs(xQ));
        sigmaPmlY = (LdomY <= abs(yQ))./(LdomY+LpmlY-abs(yQ));
        gammaPmlX = ones(size(xQ)) - sigmaPmlX/(1i*k);
        gammaPmlY = ones(size(yQ)) - sigmaPmlY/(1i*k);
        coefP = gammaPmlX.*gammaPmlY;
        coefU = gammaPmlX./gammaPmlY;
        coefV = gammaPmlY./gammaPmlX;
        matMelP = transpose(shapePhyQ) * (weightsQ .* coefP .* shapePhyQ);
        matMelU = transpose(shapePhyQ) * (weightsQ .* coefU .* shapePhyQ);
        matMelV = transpose(shapePhyQ) * (weightsQ .* coefV .* shapePhyQ);
        matIIel = [
            -1i*k*matMelP  -matDXel                          -matDYel                         ;
            -matDXel      -1i*k*matMelU                      zeros(numDofPerTRI,numDofPerTRI) ;
            -matDYel      zeros(numDofPerTRI,numDofPerTRI)  -1i*k*matMelV                     ];
    end
    
    rhsIel = [
        -1/(1i*k)*rhsPel    ;
        zeros(numDofPerTRI,1) ;
        zeros(numDofPerTRI,1) ];
    
    % ---------------------------------------------------------------------
    % Surface terms
    % ---------------------------------------------------------------------
    
    matIGel = zeros(3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);
    
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
        [xQ, yQ] = locToGloLIN(uLinQ,V1,V2);
        Jdxdu = norm(V2-V1) * 0.5;  % [ dx/du ]
        detJdxdu = abs(det(Jdxdu));
        
        % Orientation
        orientation = ones(dofm.numDofPerLIN,1);
        orientation2 = ones(dofm.numDofPerLIN,1);
        if(n1(fac) > n2(fac))
            orientation(3:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerEdg-1);
            orientation2(1:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerLIN-1);
        end
        orientation = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation);
        orientation2 = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation2);
        
        % Shape functions (physical space)
        shapePhyQ = shapePhyLinQ * orientation;
        if (BASIS == 1)
            shapeAuxQ = shapeAuxLinQ * orientation2 / sqrt(detJdxdu);
        else
            shapeAuxQ = shapePhyQ;
        end
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        weightsQ = weightsLinQ .* detJdxdu;
        
        % -----------------------------------------------------------------
        % Physical equations
        % -----------------------------------------------------------------
        
        % Local ID for interior unknowns and incoming characteristics
        idLocP = 0*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocU = 1*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocV = 2*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocG = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;
        
        % Element matrices (local element-wise system)
        matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 0.5*tau * transpose(shapePhyQ) * (weightsQ .*     shapePhyQ);
        matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + 0.5     * transpose(shapePhyQ) * (weightsQ .* nx.*shapePhyQ);
        matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + 0.5     * transpose(shapePhyQ) * (weightsQ .* ny.*shapePhyQ);
        matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 0.5     * transpose(shapePhyQ) * (weightsQ .*     shapeAuxQ);
        
        matIIel(idLocU,idLocP) = matIIel(idLocU,idLocP) + 0.5     * transpose(nx.*shapePhyQ) * (weightsQ .*     shapePhyQ);
        matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + 0.5/tau * transpose(nx.*shapePhyQ) * (weightsQ .* nx.*shapePhyQ);
        matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + 0.5/tau * transpose(nx.*shapePhyQ) * (weightsQ .* ny.*shapePhyQ);
        matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + 0.5/tau * transpose(nx.*shapePhyQ) * (weightsQ .*     shapeAuxQ);
        
        matIIel(idLocV,idLocP) = matIIel(idLocV,idLocP) + 0.5     * transpose(ny.*shapePhyQ) * (weightsQ .*     shapePhyQ);
        matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + 0.5/tau * transpose(ny.*shapePhyQ) * (weightsQ .* nx.*shapePhyQ);
        matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + 0.5/tau * transpose(ny.*shapePhyQ) * (weightsQ .* ny.*shapePhyQ);
        matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + 0.5/tau * transpose(ny.*shapePhyQ) * (weightsQ .*     shapeAuxQ);
        
        % -----------------------------------------------------------------
        % Auxiliary equations
        % -----------------------------------------------------------------
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (triNeigh > 0)
            
            % Elemental matrices (interface condition)
            matGGel = transpose(shapeAuxQ) * (weightsQ .* shapeAuxQ);
            matGIel = [ -tau*transpose(shapeAuxQ) * (weightsQ .* shapePhyQ), ...
                transpose(shapeAuxQ) * (weightsQ .* nx.*shapePhyQ), ...
                transpose(shapeAuxQ) * (weightsQ .* ny.*shapePhyQ)];
            
            % Global ID for auxiliary and exterior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            dofExt = dofm.locFacNeigh(facNeigh,:);
            idExtP = 0*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = 1*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtV = 2*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtI = [idExtP idExtU idExtV];
            
            % Assembling
            matGIx(idGloG,:) = idGloG'*ones(1,size(idExtI,2));
            matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
            matGIy(idGloG,:) = ones(size(idGloG,2),1)*idExtI;
            matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
            matGIv(idGloG,:) = matGIel;
            matGGv(idGloG,:) = matGGel;
            matGGvInv(idGloG,:) = inv(matGGel);
            
        else
            
            % Source terms
            [solQ, solDxQ, solDyQ, ~] = mySol(xQ, yQ);
            rhsPel = transpose(shapeAuxQ) * (weightsQ .* solQ);
            rhsNUel = transpose(shapeAuxQ) * (weightsQ .* (nx.*solDxQ + ny.*solDyQ)) / (1i*k);
            
            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));
            
            % Elemental matrices and RHS vectors (boundary conditions)
            matGGel = transpose(shapeAuxQ) * (weightsQ .* shapeAuxQ);
            matGIel = zeros(dofm.numDofPerLIN,3*dofm.numDofPerTRI);
            rhsGel = zeros(dofm.numDofPerLIN,1);
            switch BC
                case {'DIR0','DIR'}
                    matGIel = [tau*transpose(shapeAuxQ) * (weightsQ .* shapePhyQ), ...
                        transpose(shapeAuxQ) * (weightsQ .* nx.*shapePhyQ), ...
                        transpose(shapeAuxQ) * (weightsQ .* ny.*shapePhyQ)];
                    if(strcmp(BC,'DIR'))
                        rhsGel = +2*tau*rhsPel;
                    end
                case {'NEU0','NEU'}
                    matGIel = [-tau*transpose(shapeAuxQ) * (weightsQ .* shapePhyQ), ...
                        -transpose(shapeAuxQ) * (weightsQ .* nx.*shapePhyQ), ...
                        -transpose(shapeAuxQ) * (weightsQ .* ny.*shapePhyQ)];
                    if(strcmp(BC,'NEU'))
                        rhsGel = -2*rhsNUel;
                    end
                case {'ABC','ROB'}
                    matGIel = [tau*transpose(shapeAuxQ) * (weightsQ .* shapePhyQ), ...
                        transpose(shapeAuxQ) * (weightsQ .* nx.*shapePhyQ), ...
                        transpose(shapeAuxQ) * (weightsQ .* ny.*shapePhyQ)] * (1-tau)/(1+tau);
                    if(strcmp(BC,'ROB'))
                        rhsGel = +(rhsPel - rhsNUel) * (2*tau)/(1+tau);
                    end
                otherwise
                    error('BAD BOUNDARY CONDITION.');
            end
            
            % Global ID for auxiliary unknowns and interior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            idGloP = 0*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloU = 1*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloV = 2*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloI = [idGloP idGloU idGloV];
            
            % Assembling
            matGIx(idGloG,:) = idGloG'*ones(1,size(idGloI,2));
            matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
            matGIy(idGloG,:) = ones(size(idGloG,2),1)*idGloI;
            matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
            matGIv(idGloG,:) = matGIel;
            matGGv(idGloG,:) = matGGel;
            matGGvInv(idGloG,:) = inv(matGGel);
            rhsG(idGloG) = rhsGel;
        end
    end
    
    % ---------------------------------------------------------------------
    % Matrix assembling
    % ---------------------------------------------------------------------
    
    % Global ID of unknowns
    dofGloP = 0*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloU = 1*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloV = 2*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloI = [dofGloP dofGloU dofGloV];
    dofGloG = dofm.locToGloFAC(tri,:);
    
    % Assembling
    idTRI = (tri-1)*3*dofm.numDofPerTRI + (1:3*dofm.numDofPerTRI);
    matIIx(idTRI,:) = dofGloI'*ones(1,size(dofGloI,2));
    matIGx(idTRI,:) = dofGloI'*ones(1,size(dofGloG,2));
    matIIy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloI;
    matIGy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloG;
    matIIv(idTRI,:) = matIIel;
    matIGv(idTRI,:) = matIGel;
    matIIvInv(idTRI,:) = inv(matIIel);
    rhsI(dofGloI) = rhsIel;
    
    %condLoc(tri) = cond(full(matIIel));
    
end

% Sparse memory storage
matII    = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG    = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofFAC);
matGI    = sparse(matGIx, matGIy, matGIv, numDofFAC, 3*numDofTRI);
matGG    = sparse(matGGx, matGGy, matGGv, numDofFAC, numDofFAC);
matIIinv = sparse(matIIx, matIIy, matIIvInv, 3*numDofTRI, 3*numDofTRI);
matGGinv = sparse(matGGx, matGGy, matGGvInv, numDofFAC, numDofFAC);

% -------------------------------------------------------------------------
% Build and solve full system
% -------------------------------------------------------------------------

% Matrix partition
sysA.matII = matII;
sysA.matIG = matIG;
sysA.matGI = matGI;
sysA.matGG = matGG;
sysA.matIIinv = matIIinv;
sysA.matGGinv = matGGinv;
sysA.rhsI = rhsI;
sysA.rhsG = rhsG;

% Full system
sysA.matA = [ matII matIG ; matGI matGG ];
sysA.rhsA = [ rhsI ; rhsG ];

% Reduced system
sysA.matS = matGG - matGI*(matIIinv*matIG);
sysA.rhsS = rhsG - matGI*(matIIinv*rhsI);

% Physical system
sysA.matPhy = matII - matIG*(matGGinv*matGI);
sysA.rhsPhy = rhsI - matIG*(matGGinv*rhsG);

% Preconditionning
if (PREC == 1)
    sysA.matP = matGG;
    sysA.matPinv = matGGinv;
else
    sysA.matP = 1;
    sysA.matPinv = 1;
end

% Compute solution
solG = sysA.matS\sysA.rhsS;
solI = matIIinv*(rhsI-matIG*solG);

end