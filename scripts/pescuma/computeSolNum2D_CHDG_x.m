% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave, Simone Pescuma

function [solI, sysA] = computeSolNum2D_CHDG_x(mesh, dofm, PREC, A, B)

global k edgTagToBC
global rho c eta k

p = 0; % exponent of the power mean for the definition of \eta_F

numDofTRI = dofm.numDofTRI;
numDofFAC = dofm.numDofFAC;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);

% Shape functions and derivatives (reference space)
shapeRefLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeRefLinDuQ = functionsShapeDerLIN(uLinQ, dofm.degree);
shapeRefTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeRefTriDuQ, shapeRefTriDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% Global matrices
matIIx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIGx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIGy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIGv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIHx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIHy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIHv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIFx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIFy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIFv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);

matGIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGGx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGHx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGHy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGHv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGFx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGFy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGFv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGvInv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);

matHIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matHIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matHIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matHGx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matHGy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matHGv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matHHx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matHHy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matHHv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matHFx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matHFy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matHFv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);

matFIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matFIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matFIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matFGx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFGy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFGv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFHx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFHy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFHv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFFx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFFy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFFv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);

% Global RHS vectors
rhsI = zeros(3*numDofTRI,1);
rhsG = zeros(numDofFAC,1);
rhsH = zeros(numDofFAC,1);
rhsF = zeros(numDofFAC,1);

tic
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
    shapeLinQ = shapeRefTriQ * orientation;
    shapeTriDxQ = (shapeRefTriDuQ * Jdudx(1,1) + shapeRefTriDvQ * Jdudx(2,1)) * orientation;
    shapeTriDyQ = (shapeRefTriDuQ * Jdudx(1,2) + shapeRefTriDvQ * Jdudx(2,2)) * orientation;
    
    % Source terms
    [~, ~, ~, rhsQ, ~, ~] = mySol(xQ, yQ);
    
    % Elemental matrices and RHS vectors
    matMel = shapeLinQ' * (weightsTriQ .* shapeLinQ) * detJdxdu;
    matDXel = shapeTriDxQ' * (weightsTriQ .* shapeLinQ) * detJdxdu;
    matDYel = shapeTriDyQ' * (weightsTriQ .* shapeLinQ) * detJdxdu;
    rhsPel = shapeLinQ' * (weightsTriQ .* rhsQ) * detJdxdu;
    
    matIIel = [
        -1i*k(tri)/eta(tri)*matMel  -matDXel                          -matDYel                         ;
        -matDXel                    -1i*k(tri)*eta(tri)*matMel        zeros(numDofPerTRI,numDofPerTRI) ;
        -matDYel                    zeros(numDofPerTRI,numDofPerTRI)  -1i*k(tri)*eta(tri)*matMel       ];
    
    rhsIel = [
        -1/(1i*k(tri)*eta(tri))*rhsPel ;
        zeros(numDofPerTRI,1)          ;
        zeros(numDofPerTRI,1)          ];
    
    % ---------------------------------------------------------------------
    % Surface terms
    % ---------------------------------------------------------------------
    
    matIGel = zeros(3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);
    matIHel = zeros(3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);
    matIFel = zeros(3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);
    
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

        triNeigh = mesh.mapTriToTri(tri,fac);

        if (triNeigh>0)
            if p > 100
                etaF = max(eta(tri),eta(triNeigh));
                kF = max(k(tri),k(triNeigh));
            elseif p < -100
                etaF = min(eta(tri),eta(triNeigh));
                kF = min(k(tri),k(triNeigh));
            elseif p == 0
                etaF = sqrt(eta(tri)*eta(triNeigh));
                kF = sqrt(k(tri)*k(triNeigh));
            else
                etaF = ((eta(tri)^p+eta(triNeigh)^p)/2)^(1/p);
                kF = ((k(tri)^p+k(triNeigh)^p)/2)^(1/p);
            end
        else
            etaF = eta(tri);
            kF = k(tri);
        end

        % Mapping
        V1 = mesh.coord(n1(fac),:);
        V2 = mesh.coord(n2(fac),:);
        [xQ, yQ] = locToGloLIN(uLinQ,V1,V2);
        Jdxdu = norm(V2-V1) * 0.5;  % [ dx/du ]
        
        % Orientation
        orientation = ones(dofm.numDofPerLIN,1);
        if(n1(fac) > n2(fac))
            orientation(3:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        orientation = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation);
        
        % Shape functions (physical space)
        shapeLinQ = shapeRefLinQ * orientation;
        shapeLinDsQ = shapeRefLinDuQ * orientation;
        
        % Mass matrices (physical space)
        matM_Lin = shapeLinQ' * (weightsLinQ .* shapeLinQ) * Jdxdu;
        matK_Lin = shapeLinDsQ' * (weightsLinQ .* shapeLinDsQ) / (Jdxdu);
%         matC_Lin = 0*matM_Lin; matC_Lin(1,1) = 1; matC_Lin(2,2) = 1;
        
        if (A == 1)
            matB_Lin = matM_Lin;
        end
        if (A == 2)
            matB_Lin = matM_Lin + 0.5/(kF^2) * matK_Lin;
        end
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        % -----------------------------------------------------------------
        % Physical equations
        % -----------------------------------------------------------------
        
        % Local ID for interior unknowns and incoming characteristics
        idLocP = 0*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocU = 1*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocV = 2*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocG = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;
        idLocH = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;
        idLocF = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;
        
        % Element matrices (local element-wise system)
        matIHel(idLocP,idLocH) = matIHel(idLocP,idLocH) + 0.5/etaF * matM_Lin;
        matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 0.5/etaF * matM_Lin;
        
        matIFel(idLocU,idLocF) = matIFel(idLocU,idLocF) + nx  * matM_Lin;
        matIFel(idLocV,idLocF) = matIFel(idLocV,idLocF) + ny  * matM_Lin;
        
        % -----------------------------------------------------------------
        % Incoming characteristics
        % -----------------------------------------------------------------
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (triNeigh > 0)
            
            % Elemental matrices (interface condition)
            matGGel = matM_Lin;
            matGHel = -matM_Lin;
            
            % Global ID for auxiliary and exterior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            idLocExtH = (1:dofm.numDofPerLIN) + (facNeigh-1)*dofm.numDofPerLIN;
            idExtH = dofm.locToGloFAC(triNeigh,idLocExtH);
            idExtH([1 2]) = idExtH([2 1]);
            
            % Assembling
            matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
            matGHx(idGloG,:) = idGloG'*ones(1,size(idExtH,2));   
            matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
            matGHy(idGloG,:) = ones(size(idGloG,2),1)*idExtH;               
            matGGv(idGloG,:) = matGGel;
            matGHv(idGloG,:) = matGHel;
            matGGvInv(idGloG,:) = inv(matGGel);

        else

            % Source terms    
            [solQ, solDxQ, solDyQ, ~, ~, ~] = mySol(xQ, yQ);
            rhsPel = shapeLinQ' * (weightsLinQ .* solQ) * Jdxdu; 
            rhsUel = shapeLinQ' * (weightsLinQ .* solDxQ) * Jdxdu / (1i*kF*etaF);  
            rhsVel = shapeLinQ' * (weightsLinQ .* solDyQ) * Jdxdu / (1i*kF*etaF);

            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));
            
            % Elemental matrices and RHS vectors (boundary conditions)
            switch BC
                case 'DIR'
                    matGGel = matB_Lin;                                    % ====================================
                    matGHel = matB_Lin;                                    % ====================================
                    rhsGel  = +2*rhsPel;
                case 'NEU'
                    matGGel = matM_Lin;
                    matGHel = -matM_Lin;
                    rhsGel  = -2*etaF*(nx*rhsUel + ny*rhsVel);
                case 'ABC'
                    matGGel = matB_Lin + matM_Lin;                         % ====================================
                    matGHel = matB_Lin - matM_Lin;                         % ====================================
                    rhsGel  = +2*(rhsPel - etaF*(nx*rhsUel  + ny*rhsVel));
                case 'ROB'
                    matGGel = matB_Lin + matM_Lin;                         % ====================================
                    matGHel = matB_Lin - matM_Lin;                         % ====================================
                    rhsGel  = +2*(rhsPel - etaF*(nx*rhsUel  + ny*rhsVel));   
                otherwise
                    error('BAD BOUNDARY CONDITION.');
            end
            
            % Global ID for auxiliary unknowns and interior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            idGloH = dofm.locToGloFAC(tri,idLocH);
            
            % Assembling
            matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
            matGHx(idGloG,:) = idGloG'*ones(1,size(idGloH,2));
            matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
            matGHy(idGloG,:) = ones(size(idGloG,2),1)*idGloH;
            matGGv(idGloG,:) = matGGel;
            matGHv(idGloG,:) = matGHel;
            matGGvInv(idGloG,:) = inv(matGGel);
            
            rhsG(idGloG) = rhsGel;
            
        end

        % -----------------------------------------------------------------
        % Outgoing characteristic equations
        % -----------------------------------------------------------------
        
        % Elemental matrices (interface condition)
        if (B == 1)  % Upwind conditions (CHDG 1 & CHDG 3)
%             matHHel = matM_Lin;
%             matHIel = 2 * [ -matM_Lin, -nx*etaF*matM_Lin, -ny*etaF*matM_Lin ];
%             matHGel = - matM_Lin;
%             matHFel = 2 * matM_Lin;
            matHHel = matM_Lin;
            matHIel = [ -matM_Lin, -nx*etaF*matM_Lin, -ny*etaF*matM_Lin ];
            matHGel = 0 * matM_Lin;
            matHFel = 0 * matM_Lin;
        end
        if (B == 2)  % High-order conditions (CHDG 2)
            matHHel = matB_Lin;
            matHIel = [ -matM_Lin, -nx*etaF*matB_Lin, -ny*etaF*matB_Lin ];
            matHGel = 0 * matM_Lin;
            matHFel = 0 * matM_Lin;
        end
        
        % Global ID for flux and exterior unknowns
        idGloH = dofm.locToGloFAC(tri,idLocG);
        idGloF = dofm.locToGloFAC(tri,idLocF);
        idGloP = 0*numDofTRI + dofm.locToGloTRI(tri,idLocP);
        idGloU = 1*numDofTRI + dofm.locToGloTRI(tri,idLocP);
        idGloV = 2*numDofTRI + dofm.locToGloTRI(tri,idLocP);
        idGloI = [idGloP idGloU idGloV];
        
        % Assembling
        matHIx(idGloH,:) = idGloH'*ones(1,size(idGloI,2));
        matHHx(idGloH,:) = idGloH'*ones(1,size(idGloH,2));
        matHGx(idGloH,:) = idGloH'*ones(1,size(idGloG,2));
        matHFx(idGloH,:) = idGloH'*ones(1,size(idGloF,2));
        matHIy(idGloH,:) = ones(size(idGloH,2),1)*idGloI;
        matHHy(idGloH,:) = ones(size(idGloH,2),1)*idGloH;
        matHGy(idGloH,:) = ones(size(idGloH,2),1)*idGloG;
        matHFy(idGloH,:) = ones(size(idGloH,2),1)*idGloF;
        matHIv(idGloH,:) = matHIel;
        matHHv(idGloH,:) = matHHel;
        matHGv(idGloH,:) = matHGel;
        matHFv(idGloH,:) = matHFel;
        
        % -----------------------------------------------------------------
        % Flux equations
        % -----------------------------------------------------------------
        
        % Elemental matrices (interface condition)
        matFFel = matM_Lin;
        matFHel = - 0.5 * matB_Lin;                                        % ====================================
        matFGel = - 0.5 * matB_Lin;                                        % ====================================
        
        % Global ID for flux and exterior unknowns
        idGloG = dofm.locToGloFAC(tri,idLocG);
        idGloH = dofm.locToGloFAC(tri,idLocH);
        idGloF = dofm.locToGloFAC(tri,idLocF);
        
        % Assembling
        matFGx(idGloF,:) = idGloF'*ones(1,size(idGloG,2));
        matFHx(idGloF,:) = idGloF'*ones(1,size(idGloH,2));
        matFFx(idGloF,:) = idGloF'*ones(1,size(idGloF,2));
        matFGy(idGloF,:) = ones(size(idGloF,2),1)*idGloG;
        matFHy(idGloF,:) = ones(size(idGloF,2),1)*idGloH;
        matFFy(idGloF,:) = ones(size(idGloF,2),1)*idGloF;
        matFGv(idGloF,:) = matFGel;
        matFHv(idGloF,:) = matFHel;
        matFFv(idGloF,:) = matFFel;
        
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
    dofGloH = dofm.locToGloFAC(tri,:);
    dofGloF = dofm.locToGloFAC(tri,:);
    
    % Assembling
    idTRI = (tri-1)*3*dofm.numDofPerTRI + (1:3*dofm.numDofPerTRI);
    matIIx(idTRI,:) = dofGloI'*ones(1,size(dofGloI,2));
    matIGx(idTRI,:) = dofGloI'*ones(1,size(dofGloG,2));
    matIHx(idTRI,:) = dofGloI'*ones(1,size(dofGloH,2));
    matIFx(idTRI,:) = dofGloI'*ones(1,size(dofGloF,2));
    matIIy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloI;
    matIGy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloG;
    matIHy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloH;
    matIFy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloF;
    matIIv(idTRI,:) = matIIel;
    matIGv(idTRI,:) = matIGel;
    matIHv(idTRI,:) = matIHel;
    matIFv(idTRI,:) = matIFel;
    
    rhsI(dofGloI) = rhsIel;
    
end
toc

% Sparse memory storage
matII = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofFAC);
matIH = sparse(matIHx, matIHy, matIHv, 3*numDofTRI, numDofFAC);
matIF = sparse(matIFx, matIFy, matIFv, 3*numDofTRI, numDofFAC);
matGI = sparse(numDofFAC, 3*numDofTRI);
matGG = sparse(matGGx, matGGy, matGGv, numDofFAC, numDofFAC);
matGH = sparse(matGHx, matGHy, matGHv, numDofFAC, numDofFAC);
matGF = sparse(numDofFAC, numDofFAC);
matHI = sparse(matHIx, matHIy, matHIv, numDofFAC, 3*numDofTRI);
matHG = sparse(matHGx, matHGy, matHGv, numDofFAC, numDofFAC);
matHH = sparse(matHHx, matHHy, matHHv, numDofFAC, numDofFAC);
matHF = sparse(matHFx, matHFy, matHFv, numDofFAC, numDofFAC);
matFI = sparse(numDofFAC, 3*numDofTRI);
matFG = sparse(matFGx, matFGy, matFGv, numDofFAC, numDofFAC);
matFH = sparse(matFHx, matFHy, matFHv, numDofFAC, numDofFAC);
matFF = sparse(matFFx, matFFy, matFFv, numDofFAC, numDofFAC);

matGGinv = sparse(matGGx, matGGy, matGGvInv, numDofFAC, numDofFAC);

% -------------------------------------------------------------------------
% Build and solve full system
% -------------------------------------------------------------------------

% Matrix partition
disp('--- Matrix partition ---');
tic

X = matHF - matHH * (matFH \ matFF);
Y = matHG - matHH * (matFH \ matFG);

sysA.matII = matII + matIH * (matFH \ (matFF * (X \ matHI))) - matIF * (X \ matHI);
sysA.matIG = matIG + matIH * (matFH \ (matFF * (X \ Y)) - matFH \ matFG) - matIF * (X \ Y);
sysA.matGI = matGI + matGH * (matFH \ (matFF * (X \ matHI)));
sysA.matGG = matGG + matGH * (matFH \ (matFF * (X \ Y)) - matFH \ matFG);

sysA.matGGinv = inv(sysA.matGG);
sysA.rhsI = rhsI;
sysA.rhsG = rhsG;
toc

sysA.matIIinv = sysA.matII;
for tri=1:mesh.numTri
    dofGloP = 0*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloU = 1*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloV = 2*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloI = [dofGloP dofGloU dofGloV];
    sysA.matIIinv(dofGloI,dofGloI) = inv(sysA.matII(dofGloI,dofGloI));
end

% Reduced system
disp('--- Reduced system ---');
tic
%sysA.matS = sysA.matGG - sysA.matGI*(sysA.matII\sysA.matIG);
%sysA.rhsS = sysA.rhsG - sysA.matGI*(sysA.matII\sysA.rhsI);
sysA.matS = sysA.matGG - sysA.matGI*(sysA.matIIinv*sysA.matIG);
sysA.rhsS = sysA.rhsG - sysA.matGI*(sysA.matIIinv*sysA.rhsI);
toc

% Physical system
disp('--- Physical system ---');
tic
sysA.matPhy = sysA.matII - sysA.matIG*(sysA.matGG\sysA.matGI);
sysA.rhsPhy = sysA.rhsI - sysA.matIG*(sysA.matGG\sysA.rhsG);
toc

% Preconditionning
disp('--- Preconditionning ---');
tic
sysA.matGGinv = matGGinv;
if (PREC == 1)
    sysA.matP = matGG;
    sysA.matPinv = matGGinv;
else
    sysA.matP = 1;
    sysA.matPinv = 1;
end
toc

% Compute direct solution
disp('--- Compute direct solution ---');
tic
solG = sysA.matS\sysA.rhsS;
sol = sysA.matII\(sysA.rhsI - sysA.matIG*solG);
solI = sol(1:3*numDofTRI);
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
        error('BAD BOUNDARY TAG.')
end
end