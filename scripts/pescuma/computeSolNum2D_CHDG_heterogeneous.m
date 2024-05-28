% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave, Simone Pescuma

function [solI, sysA] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, BASIS, PREC)

global omega edgTagToBC
global rho c eta k

numDofTRI = dofm.numDofTRI;
numDofFAC = dofm.numDofFAC;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);

% Shape functions and derivatives (reference space)
shapePhyLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeAuxLinQ = functionsLegendre(uLinQ, dofm.degree);
shapePhyTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
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
matPPv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matPPvInv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);

% Global RHS vectors
rhsI = zeros(3*numDofTRI,1);
rhsG = zeros(numDofFAC,1);

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
    [~, ~, ~, rhsQ, ~, ~] = mySol(xQ, yQ);

    % Elemental matrices and RHS vectors
    matMel = shapePhyQ' * (weightsTriQ .* shapePhyQ) * detJdxdu;
    matDXel = shapeDxQ' * (weightsTriQ .* shapePhyQ) * detJdxdu;
    matDYel = shapeDyQ' * (weightsTriQ .* shapePhyQ) * detJdxdu;
    rhsPel = shapePhyQ' * (weightsTriQ .* rhsQ) * detJdxdu;
    
    matIIel = [
        -1i*k(tri)/eta(tri)*matMel  -matDXel                          -matDYel                          ;
        -matDXel                    -1i*k(tri)*eta(tri)*matMel        zeros(numDofPerTRI,numDofPerTRI)  ;
        -matDYel                    zeros(numDofPerTRI,numDofPerTRI)  -1i*k(tri)*eta(tri)*matMel        ];
    
    rhsIel = [
        -1/(1i*k(tri)*eta(tri))*rhsPel  ; 
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
        matM_IIel = shapePhyQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_IGel = shapePhyQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_GIel = shapeAuxQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_GGel = shapeAuxQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;

        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);

        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);

        % -----------------------------------------------------------------
        % Physical equations
        % -----------------------------------------------------------------
        
        % Local ID for interior unknowns and incoming characteristics
        idLocP = 0*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocU = 1*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocV = 2*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocG = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;
        
        X=0;
        if (triNeigh == 0)
            triNeigh = tri;
            X=1;
        end

        % Element matrices (local element-wise system
        matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 1/(eta(tri) + eta(triNeigh))                      * matM_IIel;
        matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + eta(tri)/(eta(tri) + eta(triNeigh))               * nx * matM_IIel;
        matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + eta(tri)/(eta(tri) + eta(triNeigh))               * ny * matM_IIel;
        matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 1/(eta(tri) + eta(triNeigh))                      * matM_IGel;
        
        matIIel(idLocU,idLocP) = matIIel(idLocU,idLocP) + eta(triNeigh)/(eta(tri) + eta(triNeigh))          * nx * matM_IIel;
        matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + eta(triNeigh)*eta(tri)/(eta(tri) + eta(triNeigh)) * nx * nx * matM_IIel;
        matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + eta(triNeigh)*eta(tri)/(eta(tri) + eta(triNeigh)) * nx * ny * matM_IIel;
        matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + eta(tri)/(eta(tri) + eta(triNeigh))               * nx * matM_IGel;

        matIIel(idLocV,idLocP) = matIIel(idLocV,idLocP) + eta(triNeigh)/(eta(tri) + eta(triNeigh))          * ny * matM_IIel;
        matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + eta(triNeigh)*eta(tri)/(eta(tri) + eta(triNeigh)) * nx * ny * matM_IIel;
        matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + eta(triNeigh)*eta(tri)/(eta(tri) + eta(triNeigh)) * ny * ny * matM_IIel;
        matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + eta(tri)/(eta(tri) + eta(triNeigh))               * ny * matM_IGel;

        if (X == 1)
            triNeigh = 0;
        end
        
        % -----------------------------------------------------------------
        % Incoming characteristics
        % -----------------------------------------------------------------
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (triNeigh > 0)
            
            % Elemental matrices (interface condition)
            matGGel = matM_GGel;
            matGIel = [-matM_GIel, eta(triNeigh)*nx*matM_GIel, eta(triNeigh)*ny*matM_GIel];
            
            % Global ID for auxiliary and exterior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            dofExt = dofm.locFacNeigh(facNeigh,:);
            idExtP = 0*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = 1*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtV = 2*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtI = [idExtP idExtU idExtV];

            % Assembling
%             coefP = max(eta(tri),eta(triNeigh))/eta(tri);
            coefP=1;

            matGIx(idGloG,:) = idGloG'*ones(1,size(idExtI,2));
            matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
            matGIy(idGloG,:) = ones(size(idGloG,2),1)*idExtI;
            matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
            matGIv(idGloG,:) = matGIel;
            matGGv(idGloG,:) = matGGel;
            matPPv(idGloG,:) = matGGel*coefP;
            matGGvInv(idGloG,:) = inv(matGGel);
            matPPvInv(idGloG,:) = inv(matGGel*coefP);

        else

            % Source terms
            [solQ, solDxQ, solDyQ, ~, ~, ~] = mySol(xQ, yQ);
            rhsPel = shapeAuxQ' * (weightsLinQ .* solQ) * Jdxdu;          
            rhsUel = shapeAuxQ' * (weightsLinQ .* solDxQ) * Jdxdu / (1i*k(tri)*eta(tri));
            rhsVel = shapeAuxQ' * (weightsLinQ .* solDyQ) * Jdxdu / (1i*k(tri)*eta(tri));

            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));

            % Elemental matrices and RHS vectors (boundary conditions)
            matGGel = matM_GGel;
            switch BC
                case 'DIR'
                    matGIel = [matM_GIel, eta(tri)*nx*matM_GIel, eta(tri)*ny*matM_GIel];
                    rhsGel  = +2*rhsPel;
                case 'NEU'
                    matGIel = [-matM_GIel, -eta(tri)*nx*matM_GIel, -eta(tri)*ny*matM_GIel];
                    rhsGel  = -2*eta(tri)*(nx*rhsUel + ny*rhsVel);
                case 'ABC'
                    matGIel = [0 .* matM_GIel, 0 .* matM_GIel, 0 .* matM_GIel];
                    rhsGel  = 0*rhsPel;
                case 'ROB'
                    matGIel = [0 .* matM_GIel, 0 .* matM_GIel, 0 .* matM_GIel];
                    rhsGel  = (rhsPel - eta(tri) * (nx*rhsUel  + ny*rhsVel));    
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
            matPPv(idGloG,:) = matGGel*coefP;
            matGGvInv(idGloG,:) = inv(matGGel);
            matPPvInv(idGloG,:) = inv(matGGel*coefP);
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

end

% Sparse memory storage
matII = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofFAC);
matGI = sparse(matGIx, matGIy, matGIv, numDofFAC, 3*numDofTRI);
matGG = sparse(matGGx, matGGy, matGGv, numDofFAC, numDofFAC);
matPP = sparse(matGGx, matGGy, matPPv, numDofFAC, numDofFAC);
matIIinv = sparse(matIIx, matIIy, matIIvInv, 3*numDofTRI, 3*numDofTRI);
matGGinv = sparse(matGGx, matGGy, matGGvInv, numDofFAC, numDofFAC);
matPPinv = sparse(matGGx, matGGy, matPPvInv, numDofFAC, numDofFAC);

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
    sysA.matP = matPP;
    sysA.matPinv = matPPinv;
else
    sysA.matP = 1;
    sysA.matPinv = 1;
end

% Compute solution
solG = sysA.matS\sysA.rhsS;
solI = matIIinv*(rhsI-matIG*solG);

end