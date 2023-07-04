% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave, Simone Pescuma

function [solI, sysA] = computeSolNum2D_CHDG_new(mesh, dofm, tau, ~, PREC, order)

global k

if order(1,1) == 2
    alfa = 1;
else
    alfa = 0;
end

if order(1,2) == 2
    beta = 1;
else
    beta = 0;
end

numDofTRI = dofm.numDofTRI;
numDofFAC = dofm.numDofFAC;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);

% Shape functions and derivatives (reference space)
shapePhyLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapePhyTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeTriDuQ, shapeTriDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);
shapePhyLinDuQ = functionsShapeDerLIN(uLinQ, dofm.degree);

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
        -1i*k*matMel  -matDXel                          -matDYel                          ;
        -matDXel      -1i*k*matMel                      zeros(numDofPerTRI,numDofPerTRI)  ;
        -matDYel      zeros(numDofPerTRI,numDofPerTRI)  -1i*k*matMel                      ];
    
    rhsIel = [
        -1/(1i*k)*rhsPel      ;
        zeros(numDofPerTRI,1) ;
        zeros(numDofPerTRI,1) ];
    
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
        shapePhyQ = shapePhyLinQ * orientation;
        shapeAuxQ = shapePhyQ;
        shapeFluQ = shapePhyQ;
        shapePhyDsQ = shapePhyLinDuQ * orientation;
        shapeAuxDsQ = shapePhyDsQ;
        shapeFluDsQ = shapePhyDsQ;
        
        % Mass matrices (physical space)
        matM_IIel = shapePhyQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_IGel = shapePhyQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_IHel = shapePhyQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_IFel = shapePhyQ' * (weightsLinQ .* shapeFluQ) * Jdxdu;

        matM_GIel = shapeAuxQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_GGel = shapeAuxQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_GHel = shapeAuxQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_GFel = shapeAuxQ' * (weightsLinQ .* shapeFluQ) * Jdxdu;
        matK_GGel = shapeAuxDsQ' * (weightsLinQ .* shapeAuxDsQ) / (Jdxdu);
        matK_GHel = shapeAuxDsQ' * (weightsLinQ .* shapeAuxDsQ) / (Jdxdu);

        matM_HIel = shapeAuxQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_HGel = shapeAuxQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_HHel = shapeAuxQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_HFel = shapeAuxQ' * (weightsLinQ .* shapeFluQ) * Jdxdu;
        matK_HIel = shapeAuxDsQ' * (weightsLinQ .* shapePhyDsQ) / (Jdxdu);
        matK_HHel = shapeAuxDsQ' * (weightsLinQ .* shapeAuxDsQ) / (Jdxdu);

        matM_FIel = shapeFluQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_FGel = shapeFluQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_FHel = shapeFluQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_FFel = shapeFluQ' * (weightsLinQ .* shapeFluQ) * Jdxdu;
        matK_FHel = shapeFluDsQ' * (weightsLinQ .* shapeAuxDsQ) / (Jdxdu);
        matK_FGel = shapeFluDsQ' * (weightsLinQ .* shapeAuxDsQ) / (Jdxdu);
        
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
        matIHel(idLocP,idLocH) = matIHel(idLocP,idLocH) + 0.5 * matM_IHel;
        matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 0.5 * matM_IGel;
        
        matIFel(idLocU,idLocF) = matIFel(idLocU,idLocF) + nx  * matM_IFel;
        matIFel(idLocV,idLocF) = matIFel(idLocV,idLocF) + ny  * matM_IFel;
        
        % -----------------------------------------------------------------
        % Incoming characteristics
        % -----------------------------------------------------------------
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (triNeigh > 0)
            
            % Elemental matrices (interface condition)
            matGGel = matM_GGel;
            matGHel = -matM_GHel;
            
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
            rhsPel = shapeAuxQ' * (weightsLinQ .* solQ) * Jdxdu; 
            rhsUel = shapeAuxQ' * (weightsLinQ .* solDxQ) * Jdxdu / (1i*k);  
            rhsVel = shapeAuxQ' * (weightsLinQ .* solDyQ) * Jdxdu / (1i*k);

            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = tagToBC(mesh.tagEdg(edgGlo));
            
            % Elemental matrices and RHS vectors (boundary conditions)
            matGGel = zeros(dofm.numDofPerLIN,dofm.numDofPerLIN);
            matGHel = zeros(dofm.numDofPerLIN,dofm.numDofPerLIN);
            rhsGel = zeros(dofm.numDofPerLIN,1);
            switch BC
                case 'DIR'
                    matGGel = matM_GGel + alfa * 0.5/(k^2) * matK_GGel;
                    matGHel = matM_GHel + alfa * 0.5/(k^2) * matK_GHel;
                    rhsGel  = +2*tau*rhsPel;
                case 'NEU'
                    matGGel = matM_GGel;
                    matGHel = -matM_GHel;
                    rhsGel  = -2*(nx*rhsUel + ny*rhsVel);
                case 'ABC'
                    matGGel = 2*matM_GGel + alfa * 0.5/(k^2) * matK_GGel;
                    matGHel = alfa * 0.5/(k^2) * matK_GHel;
                    rhsGel  = +2*(rhsPel - (nx*rhsUel  + ny*rhsVel)) * (2*tau)/(1+tau);
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
        matHHel = matM_HHel + alfa * beta * 0.5/(k^2) * matK_HHel;
        matHIel = (1 + (1 - beta)) * [-matM_HIel, -nx * (matM_HIel + alfa * beta * 0.5/(k^2)*matK_HIel), -ny *(matM_HIel + alfa * beta * 0.5/(k^2)*matK_HIel)];
        matHGel = - (1 - beta) * matM_HGel;
        matHFel = 2 * (1 - beta) * matM_HFel;

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
        matFFel = matM_FFel;         
        matFHel = - 0.5*matM_FHel - alfa * 0.25/(k^2) * matK_FHel;
        matFGel = - 0.5*matM_FGel - alfa * 0.25/(k^2) * matK_FGel; 
        
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

X = matHF - matHH * (matFH \ matFF);
Y = matHG - matHH * (matFH \ matFG);

sysA.matII = matII + matIH * (matFH \ (matFF * (X \ matHI))) - matIF * (X \ matHI);
sysA.matIG = matIG + matIH * (matFH \ (matFF * (X \ Y))) - matIF * (X \ Y) - matIH * (matFH \ matFG);
sysA.matGI = matGI + matGH * (matFH \ (matFF * (X \ matHI)));
sysA.matGG = matGG + matGH * (matFH \ (matFF * (X \ Y))) - matGH * (matFH \ matFG);
sysA.matGGinv = inv(sysA.matGG);

sysA.rhsI = rhsI;
sysA.rhsG = rhsG;

% Full system
sysA.matA = [ matII matIG matIH matIF ;
              matGI matGG matGH matGF ;
              matHI matHG matHH matHF ;
              matFI matFG matFH matFF ];
sysA.rhsA = [ rhsI ; rhsG ; rhsH ; rhsF];

% spy(sysA.matA)

% Reduced system
sysA.matS = sysA.matGG - sysA.matGI*(sysA.matII\sysA.matIG);
sysA.rhsS = sysA.rhsG - sysA.matGI*(sysA.matII\sysA.rhsI);

% Physical system
sysA.matPhy = sysA.matII - sysA.matIG*(sysA.matGG\sysA.matGI);
sysA.rhsPhy = sysA.rhsI - sysA.matIG*(sysA.matGG\sysA.rhsG);

% Preconditionning
sysA.matGGinv = matGGinv;
if (PREC == 1)
    sysA.matP = matGG;
    sysA.matPinv = matGGinv;
else
    sysA.matP = 1;
    sysA.matPinv = 1;
end

% Compute direct solution
solG = sysA.matS\sysA.rhsS;
sol = sysA.matII\(sysA.rhsI - sysA.matIG*solG);
solI = sol(1:3*numDofTRI);

% solX = sysA.matA\sysA.rhsA;
% solX = solX(1:3*numDofTRI);

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