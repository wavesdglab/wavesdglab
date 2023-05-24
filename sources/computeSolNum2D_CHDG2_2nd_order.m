% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave, Simone Pescuma

function [solI, sysA] = computeSolNum2D_CHDG2_2nd_order(mesh, dofm, tau, ~, PREC)

global k

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
shapePhyLinDuQ = functionsShapeDerLIN(uLinQ, dofm.degree);                           % NEW

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

matGIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGGx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGHx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGHy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGHv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
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

% Global RHS vectors
rhsI = zeros(3*numDofTRI,1);
rhsG = zeros(numDofFAC,1);
rhsH = zeros(numDofFAC,1);

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
        shapePhyDsQ = shapePhyLinDuQ * orientation;                          % NEW
        shapeAuxDsQ = shapePhyDsQ;                                           % NEW
        
        % Mass matrices (physical space)
        matM_IIel = shapePhyQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_IGel = shapePhyQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_IHel = shapePhyQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matK_IGel = shapePhyDsQ' * (weightsLinQ .* shapeAuxDsQ) / (Jdxdu);     % NEW
        matK_IHel = shapePhyDsQ' * (weightsLinQ .* shapeAuxDsQ) / (Jdxdu);     % NEW

        matM_GIel = shapeAuxQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_GGel = shapeAuxQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_GHel = shapeAuxQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matK_GGel = shapeAuxDsQ' * (weightsLinQ .* shapeAuxDsQ) / (Jdxdu);       % NEW
        matK_GHel = shapeAuxDsQ' * (weightsLinQ .* shapeAuxDsQ) / (Jdxdu);       % NEW

        matM_HIel = shapeAuxQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_HHel = shapeAuxQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matK_HIel = shapeAuxDsQ' * (weightsLinQ .* shapePhyDsQ) / (Jdxdu);       % NEW
        matK_HHel = shapeAuxDsQ' * (weightsLinQ .* shapeAuxDsQ) / (Jdxdu);       % NEW
        
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
        
        % Element matrices (local element-wise system)
        matIHel(idLocP,idLocH) = matIHel(idLocP,idLocH) + 0.5 * matM_IHel;
        matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 0.5 * matM_IGel;

        matIHel(idLocU,idLocH) = matIHel(idLocU,idLocH) + nx * 0.5 * matM_IHel + nx * 0.25/(k^2) * matK_IHel;
        matIHel(idLocV,idLocH) = matIHel(idLocV,idLocH) + ny * 0.5 * matM_IHel + ny * 0.25/(k^2) * matK_IHel;
        matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + nx * 0.5 * matM_IGel + nx * 0.25/(k^2) * matK_IGel;
        matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + ny * 0.5 * matM_IGel + ny * 0.25/(k^2) * matK_IGel;
        
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
            matGGel = matM_GGel;
            matGHel = zeros(dofm.numDofPerLIN,dofm.numDofPerLIN);
            rhsGel = zeros(dofm.numDofPerLIN,1);
            switch BC
                case 'DIR'
                    matGGel = matGGel + 0.5/(k^2) * matK_GGel;
                    matGHel = matM_GHel + 0.5/(k^2) * matK_GHel;
                    rhsGel  = +2*tau*rhsPel;
                case 'NEU'
                    matGHel = -matM_GHel;
                    rhsGel  = -2*(nx*rhsUel + ny*rhsVel);
                case 'ABC'
                    matGGel = matM_GGel + 0.25/(k^2) * matK_GGel;
                    matGHel = 0.25/(k^2) * matK_GHel;
                    rhsGel  = +(rhsPel - (nx*rhsUel  + ny*rhsVel)) * (2*tau)/(1+tau);
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
        matHHel = matM_HHel + 0.5/(k^2) * matK_HHel;   % UPDATED 
        matHIel = [-matM_HIel, - nx * (matM_HIel + 0.5/(k^2) * matK_HIel), - ny *(matM_HIel + 0.5/(k^2) * matK_HIel)];    % UPDATED
        
        % Global ID for flux and exterior unknowns
        idGloH = dofm.locToGloFAC(tri,idLocG);
        idGloP = 0*numDofTRI + dofm.locToGloTRI(tri,idLocP);
        idGloU = 1*numDofTRI + dofm.locToGloTRI(tri,idLocP);
        idGloV = 2*numDofTRI + dofm.locToGloTRI(tri,idLocP);
        idGloI = [idGloP idGloU idGloV];
        
        % Assembling
        matHIx(idGloH,:) = idGloH'*ones(1,size(idGloI,2));
        matHHx(idGloH,:) = idGloH'*ones(1,size(idGloH,2));
        matHIy(idGloH,:) = ones(size(idGloH,2),1)*idGloI;
        matHHy(idGloH,:) = ones(size(idGloH,2),1)*idGloH;
        matHIv(idGloH,:) = matHIel;
        matHHv(idGloH,:) = matHHel;
        
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
    
    % Assembling
    idTRI = (tri-1)*3*dofm.numDofPerTRI + (1:3*dofm.numDofPerTRI);
    matIIx(idTRI,:) = dofGloI'*ones(1,size(dofGloI,2));
    matIGx(idTRI,:) = dofGloI'*ones(1,size(dofGloG,2));
    matIHx(idTRI,:) = dofGloI'*ones(1,size(dofGloH,2));
    matIIy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloI;
    matIGy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloG;
    matIHy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloH;
    matIIv(idTRI,:) = matIIel;
    matIGv(idTRI,:) = matIGel;
    matIHv(idTRI,:) = matIHel;
    
    rhsI(dofGloI) = rhsIel;
    
end

% Sparse memory storage
matII = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofFAC);
matIH = sparse(matIHx, matIHy, matIHv, 3*numDofTRI, numDofFAC);
matGI = sparse(numDofFAC, 3*numDofTRI);
matGG = sparse(matGGx, matGGy, matGGv, numDofFAC, numDofFAC);
matGH = sparse(matGHx, matGHy, matGHv, numDofFAC, numDofFAC);
matHI = sparse(matHIx, matHIy, matHIv, numDofFAC, 3*numDofTRI);
matHG = sparse(numDofFAC, numDofFAC);
matHH = sparse(matHHx, matHHy, matHHv, numDofFAC, numDofFAC);

matGGinv = sparse(matGGx, matGGy, matGGvInv, numDofFAC, numDofFAC);

% -------------------------------------------------------------------------
% Build and solve full system
% -------------------------------------------------------------------------

% Matrix partition
sysA.matII = matII - matIH*(matHH\matHI);
sysA.matIG = matIG - matIH*(matHH\matHG);
sysA.matGI = matGI - matGH*(matHH\matHI);
sysA.matGG = matGG - matGH*(matHH\matHG);
sysA.matGGinv = inv(sysA.matGG);

% % % % % % sysA.matGGinv = matGGinv;
% sysA.matIIinv = inv(sysA.matII);

% sysA.rhsI = [rhsI; rhsH; rhsF];
% sysA.rhsG = rhsG;
% sysA.rhsI = rhsI;
% sysA.rhsG = [rhsH; rhsF; rhsG];
% sysA.rhsI = [rhsI; rhsH];
% sysA.rhsG = [rhsF; rhsG];
sysA.rhsI = rhsI;
sysA.rhsG = rhsG;

% Full system
sysA.matA = [ matII matIG matIH ;
              matGI matGG matGH ;
              matHI matHG matHH] ;
sysA.rhsA = [ rhsI ; rhsG ; rhsH];

% Reduced system
% sysA.matS = D - C*(invA*B);
% sysA.rhsS = d - C*(invA*c);
% sysA.matS = matGG - [matGI matGH matGF]*([matII matIH matIF; matHI matHH matHF; matFI matFH matFF]\[matIG; matHG; matFG]);
% sysA.rhsS = rhsG - [matGI matGH matGF]*([matII matIH matIF; matHI matHH matHF; matFI matFH matFF]\[rhsI; rhsH; rhsF]);
sysA.matS = sysA.matGG - sysA.matGI*(sysA.matII\sysA.matIG);
sysA.rhsS = sysA.rhsG - sysA.matGI*(sysA.matII\sysA.rhsI);

% Physical system
% sysA.matPhy = A - B*(invD*C);
% sysA.rhsPhy = c - B*(invD*d);
% sysA.matPhy = [matII matIH matIF; matHI matHH matHF; matFI matFH matFF] - [matIG; matHG; matFG]*(matGG\[matGI matGH matGF]);
% sysA.rhsPhy = [rhsI; rhsH; rhsF] - [matIG; matHG; matFG]*(matGG\rhsG);
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

% Compute solution
solG = sysA.matS\sysA.rhsS;
% sol = A\(c - B*solG);
% sol = [matII matIH matIF; matHI matHH matHF; matFI matFH matFF]\([rhsI; rhsH; rhsF]-[matIG; matHG; matFG]*solG);
sol = sysA.matII\(sysA.rhsI - sysA.matIG*solG);
solI = sol(1:3*numDofTRI);

% solX = sysA.matA\sysA.rhsA;   % Direct solver
% solI = solX(1:3*numDofTRI);

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