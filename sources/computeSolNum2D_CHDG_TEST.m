% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solI, sysA, condLoc] = computeSolNum2D_CHDG_TEST(mesh, dofm, tau, BASIS, PREC)

global k

numDofTRI = dofm.numDofTRI;
numDofFAC = dofm.numDofFAC;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);

% Shape functions and derivatives (reference space)

% "PHYSICAL" SHAPE FUNCTIONS
shapePhyLinQ = functionsShapeLIN(uLinQ, dofm.degree);  
shapePhyLinDuQ = functionsShapeDerLIN(uLinQ, dofm.degree);                      % NEW
shapePhyTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeTriDuQ, shapeTriDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% "AUXILIARY" SHAPE FUNCTIONS
%shapeAuxLinQ = functionsLegendre(uLinQ, dofm.degree);
shapeAuxLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeAuxLinDuQ = functionsShapeDerLIN(uLinQ, dofm.degree);                      % NEW

% "FLUX" SHAPE FUNCTIONS
% shapeFluLinQ = functionsLegendre(uLinQ, dofm.degree);                           % NEW   
shapeFluLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeFluLinDuQ = functionsShapeDerLIN(uLinQ, dofm.degree);                      % NEW

% Global matrices
matIIx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIGx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIGy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIGv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIFx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);          % NEW
matIFy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);          % NEW
matIFv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);          % NEW

matGIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGGx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGFx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);            % NEW
matGFy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);            % NEW
matGFv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);            % NEW

matFIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);          % NEW
matFIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);          % NEW
matFIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);          % NEW
matFGx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);            % NEW
matFGy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);            % NEW
matFGv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);            % NEW
matFFx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);            % NEW
matFFy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);            % NEW
matFFv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);            % NEW

% Global RHS vectors
rhsI = zeros(3*numDofTRI,1);        
rhsG = zeros(numDofFAC,1);
rhsF = zeros(numDofFAC,1);                                                     % NEW

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
    matIFel = zeros(3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);                 % NEW
    
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
        shapePhyDsQ = shapePhyLinDuQ * orientation2 / Jdxdu;                   % NEW
        
        if (BASIS == 1)
            shapeAuxQ = shapeAuxLinQ * orientation2 / sqrt(Jdxdu);
            shapeAuxDsQ = shapeAuxLinDuQ * orientation2 / sqrt(Jdxdu);         % NEW
            shapeFluQ = shapeFluLinQ * orientation2 / sqrt(Jdxdu);             % NEW
            shapeFluDsQ = shapeFluLinDuQ * orientation2 / sqrt(Jdxdu);         % NEW
        else
            shapeAuxQ = shapePhyQ;
            shapeFluQ = shapePhyQ;                                             % NEW
        end
        
        % Mass and stiffness matrices (physical space)
        matM_IIel = shapePhyQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_IGel = shapePhyQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_IFel = shapePhyQ' * (weightsLinQ .* shapeFluQ) * Jdxdu;           % NEW
        %matK_IIel = shapePhyDsQ' * (weightsLinQ .* shapePhyDsQ) * Jdxdu;       % NEW (CHECK)

        matM_GIel = shapeAuxQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_GGel = shapeAuxQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_GFel = shapeAuxQ' * (weightsLinQ .* shapeFluQ) * Jdxdu;           % NEW ( = 0 )
        %matK_GIel = shapeAuxDsQ' * (weightsLinQ .* shapePhyDsQ) * Jdxdu;       % NEW (CHECK)

        matM_FIel = shapeFluQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;           % NEW 
        matM_FGel = shapeFluQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;           % NEW
        matM_FFel = shapeFluQ' * (weightsLinQ .* shapeFluQ) * Jdxdu;           % NEW
        %matK_FIel = shapeFluDsQ' * (weightsLinQ .* shapePhyDsQ) * Jdxdu;       % NEW (CHECK)
        %matK_FFel = shapeFluDsQ' * (weightsLinQ .* shapeFluDsQ) * Jdxdu;       % NEW (CHECK)

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
        idLocF = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;                          % NEW
        
        % Element matrices (local element-wise system)  "SURFACE TERMS"
        matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 0.5*tau           * matM_IIel;
        matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + 0.5     * nx      * matM_IIel;
        matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + 0.5     * ny      * matM_IIel;

        matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 0.5               * matM_IGel;

        matIFel(idLocU,idLocF) = matIFel(idLocU,idLocF) + nx                * matM_IFel;     % NEW
        matIFel(idLocV,idLocF) = matIFel(idLocV,idLocF) + ny                * matM_IFel;     % NEW
        
        % matIIel(idLocU,idLocP) = matIIel(idLocU,idLocP) + 0.5          * nx * matM_IIel;
        % matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + 0.5/tau * nx * nx * matM_IIel;
        % matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + 0.5/tau * nx * ny * matM_IIel;
        % matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + 0.5/tau      * nx * matM_IGel;
        
        % matIIel(idLocV,idLocP) = matIIel(idLocV,idLocP) + 0.5          * ny * matM_IIel;
        % matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + 0.5/tau * nx * ny * matM_IIel;
        % matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + 0.5/tau * ny * ny * matM_IIel;
        % matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + 0.5/tau      * ny * matM_IGel;

        % matIFel(idLocU,idLocP) = matIFel(idLocU,idLocP) + nx * matM_IFel;                    % NEW
        % matIFel(idLocU,idLocU) = matIFel(idLocU,idLocU) + (nx * matM_IFel).*0;               % NEW
        % matIFel(idLocU,idLocV) = matIFel(idLocU,idLocV) + (nx * matM_IFel).*0;               % NEW
  
        % matIFel(idLocV,idLocP) = matIFel(idLocU,idLocP) + ny * matM_IFel;                    % NEW
        % matIFel(idLocV,idLocU) = matIFel(idLocU,idLocU) + (ny * matM_IFel).*0;               % NEW
        % matIFel(idLocV,idLocV) = matIFel(idLocU,idLocV) + (ny * matM_IFel).*0;               % NEW

        % -----------------------------------------------------------------
        % Auxiliary equations
        % -----------------------------------------------------------------
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (triNeigh > 0)
            
            % Elemental matrices (interface condition)
            matGGel = matM_GGel;
            matGIel = [-tau*matM_GIel, nx*matM_GIel, ny*matM_GIel];           % UPDATED
            matGFel = matM_GFel;                                              % NEW

            % Global ID for auxiliary and exterior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            idGloF = dofm.locToGloFAC(tri,idLocF);                            % NEW
            dofExt = dofm.locFacNeigh(facNeigh,:);
            idExtP = 0*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = 1*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtV = 2*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtI = [idExtP idExtU idExtV];
            
            % Assembling
            matGIx(idGloG,:) = idGloG'*ones(1,size(idExtI,2));
            matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
            matGFx(idGloG,:) = idGloG'*ones(1,size(idGloF,2));                 % NEW
            matGIy(idGloG,:) = ones(size(idGloG,2),1)*idExtI;
            matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
            matGFy(idGloG,:) = ones(size(idGloG,2),1)*idGloF;
            matGIv(idGloG,:) = matGIel;
            matGGv(idGloG,:) = matGGel;
            matGFv(idGloG,:) = 0.*matGFel;                                     % NEW

            matGGvInv(idGloG,:) = inv(matGGel);

        else
            
            % Source terms
            [solQ, solDxQ, solDyQ, ~] = mySol(xQ, yQ);
            rhsPel = shapeAuxQ' * (weightsLinQ .* solQ) * Jdxdu;
            rhsUel = shapeAuxQ' * (weightsLinQ .* solDxQ) * Jdxdu / (1i*k);
            rhsVel = shapeAuxQ' * (weightsLinQ .* solDyQ) * Jdxdu / (1i*k);
            
            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = tagToBC(mesh.tagEdg(edgGlo));
            
            % Elemental matrices and RHS vectors (boundary conditions)    
            matGGel = matM_GGel;
            matGIel = zeros(dofm.numDofPerLIN,3*dofm.numDofPerTRI);
            matGFel = 0.*matM_GFel;                                          % NEW
            rhsGel = zeros(dofm.numDofPerLIN,1);
            switch BC
                case 'DIR'
                    matGIel = [+tau*matM_GIel, +nx*matM_GIel, +ny*matM_GIel];    
                    rhsGel  = +2*tau*rhsPel;
                case 'NEU'
                    matGIel = [-tau*matM_GIel, -nx*matM_GIel, -ny*matM_GIel];   
                    rhsGel  = -2*(nx*rhsUel + ny*rhsVel);
                case 'ABC'
                    matGIel = [+tau*matM_GIel, +nx*matM_GIel, +ny*matM_GIel] * (1-tau)/(1+tau);  % = 0
                    rhsGel  = +(rhsPel - (nx*rhsUel + ny*rhsVel)) * (2*tau)/(1+tau);
                otherwise
                    error('BAD BOUNDARY CONDITION.');
            end
            
            % Global ID for auxiliary unknowns and interior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            idGloF = dofm.locToGloFAC(tri,idLocF);                            % NEW
            idGloP = 0*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloU = 1*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloV = 2*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloI = [idGloP idGloU idGloV]; 

            % Assembling
            matGIx(idGloG,:) = idGloG'*ones(1,size(idGloI,2));
            matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
            matGFx(idGloG,:) = idGloG'*ones(1,size(idGloF,2));                % NEW
            matGIy(idGloG,:) = ones(size(idGloG,2),1)*idGloI;
            matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
            matGFy(idGloG,:) = ones(size(idGloG,2),1)*idGloF;                 % NEW
            matGIv(idGloG,:) = matGIel;
            matGGv(idGloG,:) = matGGel;
            matGFv(idGloG,:) = matGFel;                                       % NEW

            matGGvInv(idGloG,:) = inv(matGGel);

            rhsG(idGloG) = rhsGel;

        end
    
        % -----------------------------------------------------------------
        % Flux equations
        % -----------------------------------------------------------------

        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (triNeigh > 0)
            
            % Elemental matrices (interface condition)
            matFFel = matM_FFel;                                                    % NEW
            matFIel = [-0.5*matM_FIel, -0.5*nx*matM_FIel, -0.5*ny*matM_FIel];       % NEW
            matFGel = -0.5*matM_FGel;                                               % NEW

            % Global ID for flux and exterior unknowns 
            idGloG = dofm.locToGloFAC(tri,idLocG);
            idGloF = dofm.locToGloFAC(tri,idLocF);                                  % NEW
            dofExt = dofm.locFacNeigh(facNeigh,:);
            idExtP = 0*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = 1*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtV = 2*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtI = [idExtP idExtU idExtV];
            
            % Assembling                               
            matFIx(idGloF,:) = idGloF'*ones(1,size(idExtI,2));
            matFGx(idGloF,:) = idGloF'*ones(1,size(idGloG,2));
            matFFx(idGloF,:) = idGloF'*ones(1,size(idGloF,2));
            matFIy(idGloF,:) = ones(size(idGloF,2),1)*idExtI;
            matFGy(idGloF,:) = ones(size(idGloF,2),1)*idGloG;
            matFFy(idGloF,:) = ones(size(idGloF,2),1)*idGloF;
            matFIv(idGloF,:) = matFIel;
            matFGv(idGloF,:) = matFGel;
            matFFv(idGloF,:) = matFFel;
            
            matFGvInv(idGloF,:) = inv(matFGel);

        else

            % Source terms
            [solQ, solDxQ, solDyQ, ~] = mySol(xQ, yQ);
            rhsPel = (shapeFluQ' * (weightsLinQ .* solQ) * Jdxdu);             % NEW
            rhsUel = (shapeFluQ' * (weightsLinQ .* solDxQ) * Jdxdu);           % NEW
            rhsVel = (shapeFluQ' * (weightsLinQ .* solDyQ) * Jdxdu);           % NEW
            
            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = tagToBC(mesh.tagEdg(edgGlo));
            
            % Elemental matrices and RHS vectors (boundary conditions) 
            % FIX B.C. 
             
            matFGel = -0.5*matM_FGel;                                          % NEW
            matFIel = zeros(dofm.numDofPerLIN,3*dofm.numDofPerTRI);            % NEW
            matFFel = 0*matM_FFel;                                             % NEW
            rhsFel = zeros(dofm.numDofPerLIN,1);                               % NEW
            switch BC
                case 'DIR'
                    matFIel = [-0.5*matM_FIel, -0.5*nx*matM_FIel, -0.5*ny*matM_FIel];
                    rhsFel = - rhsPel;
                    %rhsFel = +2*tau*rhsPel; % same as above, but *(-2)
                case 'NEU'
                    matFIel = [+0.5*matM_FIel, +0.5*nx*matM_FIel, +0.5*ny*matM_FIel];
                    rhsFel = nx*rhsUel + ny*rhsVel;
                    %rhsFel = -2*(nx*rhsUel + ny*rhsVel); % same as above, but *(-2)
                case 'ABC'
                    matFIel = [-0.5*matM_FIel, -0.5*nx*matM_FIel, -0.5*ny*matM_FIel] * (1-tau)/(1+tau); % = 0
                    rhsFel = -1/2 * (rhsPel - (nx*rhsUel + ny*rhsVel)) * (2*tau)/(1+tau); 
                    %rhsFel = (rhsPel - (nx*rhsUel + ny*rhsVel)) * (2*tau)/(1+tau); % same as above, but *(-2)
                otherwise
                    error('BAD BOUNDARY CONDITION.');
            end
            
            % Global ID for flux unknowns and interior unknowns
            idGloF = dofm.locToGloFAC(tri,idLocF);
            idGloP = 0*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloU = 1*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloV = 2*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloI = [idGloP idGloU idGloV];
            
            % Assembling
            matFIx(idGloF,:) = idGloF'*ones(1,size(idGloI,2));                % NEW
            matFGx(idGloF,:) = idGloF'*ones(1,size(idGloG,2));                % NEW
            matFFx(idGloF,:) = idGloF'*ones(1,size(idGloF,2));                % NEW
            matFIy(idGloF,:) = ones(size(idGloF,2),1)*idGloI;                 % NEW
            matFGy(idGloF,:) = ones(size(idGloF,2),1)*idGloG;                 % NEW
            matFFy(idGloF,:) = ones(size(idGloF,2),1)*idGloF;                 % NEW
            matFIv(idGloF,:) = matFIel;                                       % NEW
            matFGv(idGloF,:) = matFGel;                                       % NEW
            matFFv(idGloF,:) = matFFel;                                       % NEW

            matFGvInv(idGloF,:) = inv(matFGel);                               % NEW

            rhsF(idGloF) = rhsFel;

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
    dofGloF = dofm.locToGloFAC(tri,:);                                        % NEW
    
    % Assembling
    idTRI = (tri-1)*3*dofm.numDofPerTRI + (1:3*dofm.numDofPerTRI);
    matIIx(idTRI,:) = dofGloI'*ones(1,size(dofGloI,2));
    matIGx(idTRI,:) = dofGloI'*ones(1,size(dofGloG,2));
    matIFx(idTRI,:) = dofGloI'*ones(1,size(dofGloF,2));                       % NEW
    matIIy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloI;                               
    matIGy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloG;
    matIFy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloF;                        % NEW
    matIIv(idTRI,:) = matIIel;
    matIGv(idTRI,:) = matIGel;
    matIFv(idTRI,:) = matIFel;                                                % NEW

    matIIvInv(idTRI,:) = inv(matIIel);

    rhsI(dofGloI) = rhsIel;
    
    condLoc(tri) = cond(full(matIIel));
    
end

% Sparse memory storage
matII    = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG    = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofFAC);
matIF    = sparse(matIFx, matIFy, matIFv, 3*numDofTRI, numDofFAC);           % NEW
matGI    = sparse(matGIx, matGIy, matGIv, numDofFAC, 3*numDofTRI);
matGG    = sparse(matGGx, matGGy, matGGv, numDofFAC, numDofFAC);
matGF    = 0.*sparse(matGFx, matGFy, matGFv, numDofFAC, numDofFAC);          % NEW
matFI    = sparse(matFIx, matFIy, matFIv, numDofFAC, 3*numDofTRI);           % NEW
matFG    = sparse(matFGx, matFGy, matFGv, numDofFAC, numDofFAC);             % NEW
matFF    = sparse(matFFx, matFFy, matFFv, numDofFAC, numDofFAC);             % NEW

matIIinv = sparse(matIIx, matIIy, matIIvInv, 3*numDofTRI, 3*numDofTRI);
matGGinv = sparse(matGGx, matGGy, matGGvInv, numDofFAC, numDofFAC);

% -------------------------------------------------------------------------
% Build and solve full system
% -------------------------------------------------------------------------

% Matrix partition
sysA.matII = matII;
sysA.matIG = matIG;
sysA.matIF = matIF;
sysA.matGI = matGI;
sysA.matGG = matGG;
sysA.matGF = matGF;         % NEW 
sysA.matFI = matFI;         % NEW 
sysA.matFG = matFG;         % NEW 
sysA.matFF = matFF;         % NEW 
sysA.matIIinv = matIIinv;
sysA.matGGinv = matGGinv;
sysA.rhsI = rhsI;
sysA.rhsG = rhsG;
sysA.rhsF = rhsF;           % NEW 

% Full system
sysA.matA = [ matII matIG matIF; matGI matGG matGF; matFI matFG matFF ];       % UPDATED
sysA.rhsA = [ rhsI ; rhsG ; rhsF];                                             % UPDATED

% Reduced system
sysA.matS = matGG - matGI*(matIIinv*matIG);
sysA.rhsS = rhsG - matGI*(matIIinv*rhsI);

% Physical system
sysA.matPhy = matII - matIG*(matGGinv*matGI);
sysA.rhsPhy = rhsI - matIG*(matGGinv*rhsG);

spy(sysA.matA)

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
% solI = matIIinv*(rhsI-matIG*solG);
solI = sysA.matA\sysA.rhsA;

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