% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solI, sysA, condLoc] = computeSolNum2D_HDG(mesh, dofm, tau, BASIS, PREC)

global k edgTagToBC
global LdomX LdomY LpmlX LpmlY Rdom Rpml

numDofTRI = dofm.numDofTRI;
numDofLIN = dofm.numDofLIN;
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
matIGx = zeros(3*dofm.numDofPerTRI, mesh.numTri*3*dofm.numDofPerLIN);
matIGy = zeros(3*dofm.numDofPerTRI, mesh.numTri*3*dofm.numDofPerLIN);
matIGv = zeros(3*dofm.numDofPerTRI, mesh.numTri*3*dofm.numDofPerLIN);
matGIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGGx = zeros(numDofLIN,dofm.numDofPerLIN);
matGGy = zeros(numDofLIN,dofm.numDofPerLIN);
matGGv = zeros(numDofLIN,dofm.numDofPerLIN);
matIIvInv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);

% Global RHS vectors
rhsI = zeros(3*numDofTRI,1);
rhsG = zeros(numDofLIN,1);

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
    
    % Source terms
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
    
    % PML stretching (rectangular PML)
    if(~isempty(LdomX) && ~isempty(LdomY))
        if ((mean(abs(xQ)) >= LdomX) || (mean(abs(yQ)) >= LdomY))
            sigmaPmlX = (LdomX <= abs(xQ))./(LdomX+LpmlX-abs(xQ));
            sigmaPmlY = (LdomY <= abs(yQ))./(LdomY+LpmlY-abs(yQ));
            gammaPmlX = ones(size(xQ)) - sigmaPmlX/(1i*k);
            gammaPmlY = ones(size(yQ)) - sigmaPmlY/(1i*k);
            coefPml = gammaPmlX.*gammaPmlY;
            tensPmlInvXX = gammaPmlX./gammaPmlY;
            tensPmlInvYY = gammaPmlY./gammaPmlX;
            matMelP = transpose(shapePhyQ) * (weightsQ .* coefPml .* shapePhyQ);
            matMelU = transpose(shapePhyQ) * (weightsQ .* tensPmlInvXX .* shapePhyQ);
            matMelV = transpose(shapePhyQ) * (weightsQ .* tensPmlInvYY .* shapePhyQ);
            matIIel = [
                -1i*k*matMelP  -matDXel                          -matDYel                         ;
                -matDXel      -1i*k*matMelU                      zeros(numDofPerTRI,numDofPerTRI) ;
                -matDYel      zeros(numDofPerTRI,numDofPerTRI)  -1i*k*matMelV                     ];
        end
    end
    
    % PML stretching (circular PML)
    if(~isempty(Rdom))
        rQ = sqrt(xQ.*xQ + yQ.*yQ);
        if (mean(rQ) >= Rdom)
            cosT = xQ./rQ;
            sinT = yQ./rQ;
            sigmaPml = 1./(Rpml-(rQ-Rdom));
            sigmaPmlInt = -log(1-(rQ-Rdom)/Rpml);
            gammaPmlR = ones(size(rQ)) - sigmaPml/(1i*k);
            gammaPmlT = ones(size(rQ)) - sigmaPmlInt/(1i*k)./rQ;
            coefPml = gammaPmlR.*gammaPmlT;
            tensPmlInvXX = (gammaPmlR./gammaPmlT) .* cosT.*cosT + (gammaPmlT./gammaPmlR) .* (sinT.*sinT);
            tensPmlInvXY = (gammaPmlR./gammaPmlT) .* cosT.*sinT - (gammaPmlT./gammaPmlR) .* (cosT.*sinT);
            tensPmlInvYY = (gammaPmlR./gammaPmlT) .* sinT.*sinT + (gammaPmlT./gammaPmlR) .* (cosT.*cosT);
            matMelP = transpose(shapePhyQ) * (weightsQ .* coefPml .* shapePhyQ);
            matMelXX = transpose(shapePhyQ) * (weightsQ .* tensPmlInvXX .* shapePhyQ);
            matMelXY = transpose(shapePhyQ) * (weightsQ .* tensPmlInvXY .* shapePhyQ);
            matMelYY = transpose(shapePhyQ) * (weightsQ .* tensPmlInvYY .* shapePhyQ);
            matIIel = [
                -1i*k*matMelP  -matDXel        -matDYel       ;
                -matDXel       -1i*k*matMelXX  -1i*k*matMelXY ;
                -matDYel       -1i*k*matMelXY  -1i*k*matMelYY ];
        end
    end
    
    rhsIel = [
        -1/(1i*k)*rhsPel    ;
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
        
        % Mapping
        V1 = mesh.coord(n1(fac),:);
        V2 = mesh.coord(n2(fac),:);
        [xQ, yQ] = locToGloLIN(uLinQ, V1, V2);
        Jdxdu = norm(V2-V1) * 0.5;  % [ dx/du ]
        
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
            shapeAuxQ = shapeAuxLinQ * orientation2 / sqrt(Jdxdu);
        else
            shapeAuxQ = shapePhyQ;
        end
        
        % Mass matrices (physical space)
        matM_IIel = transpose(shapePhyQ) * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_IGel = transpose(shapePhyQ) * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_GIel = transpose(shapeAuxQ) * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_GGel = transpose(shapeAuxQ) * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        % -----------------------------------------------------------------
        % Physical equations
        % -----------------------------------------------------------------
        
        % Local ID for interior unknowns
        dofLocP = 0*numDofPerTRI + dofm.locFac(fac,:);
        dofLocU = 1*numDofPerTRI + dofm.locFac(fac,:);
        dofLocV = 2*numDofPerTRI + dofm.locFac(fac,:);
        dofLocI = [dofLocP, dofLocU, dofLocV];
        
        % Element matrices (local element-wise system)
        matIIel(dofLocP,dofLocP) = matIIel(dofLocP,dofLocP) + tau*matM_IIel;
        matIIel(dofLocP,dofLocU) = matIIel(dofLocP,dofLocU) + nx*matM_IIel;
        matIIel(dofLocP,dofLocV) = matIIel(dofLocP,dofLocV) + ny*matM_IIel;
        matIGel = zeros(3*dofm.numDofPerTRI,dofm.numDofPerLIN);
        matIGel(dofLocP,:) = -tau * matM_IGel;
        matIGel(dofLocU,:) =  nx  * matM_IGel;
        matIGel(dofLocV,:) =  ny  * matM_IGel;
        
        % -----------------------------------------------------------------
        % Auxiliary equations
        % -----------------------------------------------------------------
        
        matGGel = zeros(dofm.numDofPerLIN, dofm.numDofPerLIN);
        matGIel = zeros(dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
        rhsGel  = zeros(dofm.numDofPerLIN, 1);
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        
        if (triNeigh > 0)
            
            % Elemental matrices (interface condition)
            matGGel = 0.5 * matM_GGel;
            matGIel = -0.5/tau * [tau*matM_GIel, nx*matM_GIel, ny*matM_GIel];
            
        else
            
            % Source terms
            [solQ, solDxQ, solDyQ, ~] = mySol(xQ, yQ);
            rhsPel = transpose(shapeAuxQ) * (weightsLinQ .* solQ) * Jdxdu;
            rhsUel = transpose(shapeAuxQ) * (weightsLinQ .* solDxQ) * Jdxdu / (1i*k);
            rhsVel = transpose(shapeAuxQ) * (weightsLinQ .* solDyQ) * Jdxdu / (1i*k);
            
            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));
            
            % Elemental matrices and RHS vectors (boundary conditions)
            matGGel = matM_GGel;
            switch BC
                case 'DIR0'
                case 'DIR'
                    rhsGel = rhsPel;
                case 'NEU0'
                    matGIel = -1/tau * [tau*matM_GIel, nx*matM_GIel, ny*matM_GIel];
                case 'NEU'
                    matGIel = -1/tau * [tau*matM_GIel, nx*matM_GIel, ny*matM_GIel];
                    rhsGel = -(nx*rhsUel + ny*rhsVel) / tau;
                case 'ABC'
                    matGIel = -1/(1+tau) * [tau*matM_GIel, nx*matM_GIel, ny*matM_GIel];
                case 'ROB'
                    matGIel = -1/(1+tau) * [tau*matM_GIel, nx*matM_GIel, ny*matM_GIel];
                    rhsGel = (rhsPel - nx*rhsUel - ny*rhsVel) / (1+tau);
                otherwise
                    error('BAD BOUNDARY CONDITION.');
            end
        end
        
        % Global ID for edge unknowns
        edgGlo = abs(mesh.mapTriToEdg(tri,fac));
        dofGloG = dofm.locToGloLIN(edgGlo,:);
        dofLocG = 1:dofm.numDofPerLIN;
        if (BASIS ~= 1)
            if(mesh.mapTriToEdg(tri,fac) < 0)
                dofLocG(1) = 2;
                dofLocG(2) = 1;
            end
        end
        
        % -----------------------------------------------------------------
        % Matrix assembling
        % -----------------------------------------------------------------
        
        idLIN = (tri-1)*3*dofm.numDofPerLIN + (fac-1)*dofm.numDofPerLIN + (1:dofm.numDofPerLIN);
        matIGx(:,idLIN) = dofGloI'*ones(1,size(dofGloG,2));
        matIGy(:,idLIN) = ones(size(dofGloI,2),1)*dofGloG;
        matIGv(:,idLIN) = matIGel(:,dofLocG);
        matGIx(idLIN,:) = dofGloG'*ones(1,size(dofGloI(dofLocI),2));
        matGIy(idLIN,:) = ones(size(dofGloG,2),1)*dofGloI(dofLocI);
        matGIv(idLIN,:) = matGIel(dofLocG,:);
        matGGx(dofGloG,:) = dofGloG'*ones(1,size(dofGloG,2));
        matGGy(dofGloG,:) = ones(size(dofGloG,2),1)*dofGloG;
        matGGv(dofGloG,:) = matGGv(dofGloG,:) + matGGel(dofLocG,dofLocG);
        rhsG(dofGloG) = rhsG(dofGloG) + rhsGel(dofLocG);
    end
    
    % -------------------------------------------------------------------------
    % Matrix assembling
    % -------------------------------------------------------------------------
    
    idTRI = (tri-1)*3*dofm.numDofPerTRI + (1:3*dofm.numDofPerTRI);
    matIIx(idTRI,:) = dofGloI'*ones(1,size(dofGloI,2));
    matIIy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloI;
    matIIv(idTRI,:) = matIIel;
    matIIvInv(idTRI,:) = inv(matIIel);
    rhsI(dofGloI) = rhsIel;
    
    condLoc(tri) = cond(full(matIIel));
    
end

matII    = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG    = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofLIN);
matGI    = sparse(matGIx, matGIy, matGIv, numDofLIN, 3*numDofTRI);
matGG    = sparse(matGGx, matGGy, matGGv, numDofLIN, numDofLIN);
matIIinv = sparse(matIIx, matIIy, matIIvInv, 3*numDofTRI, 3*numDofTRI);

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

% Matrix partition
sysA.matII = matII;
sysA.matIG = matIG;
sysA.matGI = matGI;
sysA.matGG = matGG;
sysA.matIIinv = matIIinv;
sysA.matGGinv = inv(matGG);
sysA.rhsI = rhsI;
sysA.rhsG = rhsG;

% Full system
sysA.matA = [ matII matIG ; matGI matGG ];
sysA.rhsA = [ rhsI ; rhsG ];

% Reduced system
sysA.matS = matGG - matGI*(matIIinv*matIG);
sysA.rhsS = rhsG - matGI*(matIIinv*rhsI);

% Physical system
sysA.matPhy = matII - matIG*(sysA.matGGinv*matGI);
sysA.rhsPhy = rhsI - matIG*(sysA.matGGinv*rhsG);

% Preconditionning
if (PREC == 1)
    sysA.matP = matGG;
    sysA.matPinv = sysA.matGGinv;
else
    sysA.matP = 1;
    sysA.matPinv = 1;
end

% Compute solution
solG = sysA.matS\sysA.rhsS;
solI = matIIinv*(rhsI-matIG*solG);

end