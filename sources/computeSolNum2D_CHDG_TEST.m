% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave, Simone Pescuma

function [solI, sysA] = computeSolNum2D_CHDG_TEST(mesh, dofm, tau, ~, ~)

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

% Global matrices
matIIx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIGx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIGy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIGv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIFx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIFy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIFv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);

matGIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGGx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGFx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGFy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGFv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);

matFIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matFIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matFIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matFGx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFGy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFGv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFFx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFFy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matFFv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);

% Global RHS vectors
rhsI = zeros(3*numDofTRI,1);    
rhsG = zeros(numDofFAC,1);
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
    Jdudx = inv(Jdxdu);   % [ du/dx du/dy ; dv/dx dv/dy ]
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
        
        % Mass matrices (physical space)
        matM_IIel = shapePhyQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_IGel = shapePhyQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_IFel = shapePhyQ' * (weightsLinQ .* shapeFluQ) * Jdxdu;
        matM_GIel = shapeAuxQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_GGel = shapeAuxQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_GFel = shapeAuxQ' * (weightsLinQ .* shapeFluQ) * Jdxdu;
        matM_FIel = shapeFluQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;  
        matM_FGel = shapeFluQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu; 
        matM_FFel = shapeFluQ' * (weightsLinQ .* shapeFluQ) * Jdxdu; 

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
        idLocF = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;                
        
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

        % matIFel(idLocU,idLocP) = matIFel(idLocU,idLocP) + nx * matM_IFel;          
        % matIFel(idLocU,idLocU) = matIFel(idLocU,idLocU) + (nx * matM_IFel).*0;     
        % matIFel(idLocU,idLocV) = matIFel(idLocU,idLocV) + (nx * matM_IFel).*0;     
  
        % matIFel(idLocV,idLocP) = matIFel(idLocU,idLocP) + ny * matM_IFel;          
        % matIFel(idLocV,idLocU) = matIFel(idLocU,idLocU) + (ny * matM_IFel).*0;     
        % matIFel(idLocV,idLocV) = matIFel(idLocU,idLocV) + (ny * matM_IFel).*0;     

        % -----------------------------------------------------------------
        % Auxiliary equations
        % -----------------------------------------------------------------
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (triNeigh > 0)
            
            % Elemental matrices (interface condition)
            matGGel = matM_GGel;
            matGIel = [-tau*matM_GIel, nx*matM_GIel, ny*matM_GIel];
            matGFel = 0.*matM_GFel;

            % Global ID for auxiliary and exterior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            idGloF = dofm.locToGloFAC(tri,idLocF);
            dofExt = dofm.locFacNeigh(facNeigh,:);
            idExtP = 0*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = 1*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtV = 2*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtI = [idExtP idExtU idExtV];
            
            % Assembling
            matGIx(idGloG,:) = idGloG'*ones(1,size(idExtI,2));
            matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
            matGFx(idGloG,:) = idGloG'*ones(1,size(idGloF,2));     
            matGIy(idGloG,:) = ones(size(idGloG,2),1)*idExtI;
            matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
            matGFy(idGloG,:) = ones(size(idGloG,2),1)*idGloF;
            matGIv(idGloG,:) = matGIel;
            matGGv(idGloG,:) = matGGel;
            matGFv(idGloG,:) = matGFel;

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
            matGFel = 0.*matM_GFel; 
            rhsGel = zeros(dofm.numDofPerLIN,1);
            switch BC
                case 'DIR'
                    matGIel = [+tau*matM_GIel, +nx*matM_GIel, +ny*matM_GIel];  
                    rhsGel  = +2*tau*rhsPel;
                case 'NEU'
                    matGIel = [-tau*matM_GIel, -nx*matM_GIel, -ny*matM_GIel]; 
                    rhsGel  = -2*(nx*rhsUel + ny*rhsVel);
                case 'ABC'
                    matGIel = [+tau*matM_GIel, +nx*matM_GIel, +ny*matM_GIel] * (1-tau)/(1+tau);
                    rhsGel  = +(rhsPel - (nx*rhsUel + ny*rhsVel)) * (2*tau)/(1+tau);
                otherwise
                    error('BAD BOUNDARY CONDITION.');
            end
            
            % Global ID for auxiliary unknowns and interior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            idGloF = dofm.locToGloFAC(tri,idLocF);  
            idGloP = 0*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloU = 1*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloV = 2*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloI = [idGloP idGloU idGloV]; 

            % Assembling
            matGIx(idGloG,:) = idGloG'*ones(1,size(idGloI,2));
            matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
            matGFx(idGloG,:) = idGloG'*ones(1,size(idGloF,2));    
            matGIy(idGloG,:) = ones(size(idGloG,2),1)*idGloI;
            matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
            matGFy(idGloG,:) = ones(size(idGloG,2),1)*idGloF;       
            matGIv(idGloG,:) = matGIel;
            matGGv(idGloG,:) = matGGel;
            matGFv(idGloG,:) = matGFel;

            rhsG(idGloG) = rhsGel;

        end
    
        % -----------------------------------------------------------------
        % Flux equations
        % -----------------------------------------------------------------
        
        % Elemental matrices (interface condition)
        matFFel = matM_FFel;
        matFIel = [-0.5*matM_FIel, -0.5*nx*matM_FIel, -0.5*ny*matM_FIel];
        matFGel = -0.5*matM_FGel;
        
        % Global ID for flux and exterior unknowns
        idGloG = dofm.locToGloFAC(tri,idLocG);
        idGloF = dofm.locToGloFAC(tri,idLocF);
        idGloP = 0*numDofTRI + dofm.locToGloTRI(tri,idLocP);
        idGloU = 1*numDofTRI + dofm.locToGloTRI(tri,idLocP);
        idGloV = 2*numDofTRI + dofm.locToGloTRI(tri,idLocP);
        idGloI = [idGloP idGloU idGloV];
        
        % Assembling
        matFIx(idGloF,:) = idGloF'*ones(1,size(idGloI,2));
        matFGx(idGloF,:) = idGloF'*ones(1,size(idGloG,2));
        matFFx(idGloF,:) = idGloF'*ones(1,size(idGloF,2));
        matFIy(idGloF,:) = ones(size(idGloF,2),1)*idGloI;
        matFGy(idGloF,:) = ones(size(idGloF,2),1)*idGloG;
        matFFy(idGloF,:) = ones(size(idGloF,2),1)*idGloF;
        matFIv(idGloF,:) = matFIel;
        matFGv(idGloF,:) = matFGel;
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
    dofGloF = dofm.locToGloFAC(tri,:);
    
    % Assembling
    idTRI = (tri-1)*3*dofm.numDofPerTRI + (1:3*dofm.numDofPerTRI);
    matIIx(idTRI,:) = dofGloI'*ones(1,size(dofGloI,2));
    matIGx(idTRI,:) = dofGloI'*ones(1,size(dofGloG,2));
    matIFx(idTRI,:) = dofGloI'*ones(1,size(dofGloF,2));    
    matIIy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloI;
    matIGy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloG;
    matIFy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloF;              
    matIIv(idTRI,:) = matIIel;
    matIGv(idTRI,:) = matIGel;
    matIFv(idTRI,:) = matIFel;

    rhsI(dofGloI) = rhsIel;
    
end

% Sparse memory storage
matII    = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG    = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofFAC);
matIF    = sparse(matIFx, matIFy, matIFv, 3*numDofTRI, numDofFAC); 
matGI    = sparse(matGIx, matGIy, matGIv, numDofFAC, 3*numDofTRI);
matGG    = sparse(matGGx, matGGy, matGGv, numDofFAC, numDofFAC);
matGF    = 0.*sparse(matGFx, matGFy, matGFv, numDofFAC, numDofFAC);
matFI    = sparse(matFIx, matFIy, matFIv, numDofFAC, 3*numDofTRI); 
matFG    = sparse(matFGx, matFGy, matFGv, numDofFAC, numDofFAC); 
matFF    = sparse(matFFx, matFFy, matFFv, numDofFAC, numDofFAC);

% -------------------------------------------------------------------------
% Build and solve full system
% -------------------------------------------------------------------------

% Matrix partition
sysA.matII = matII;
sysA.matIG = matIG;
sysA.matIF = matIF;
sysA.matGI = matGI;
sysA.matGG = matGG;
sysA.matGF = matGF;
sysA.matFI = matFI;
sysA.matFG = matFG;
sysA.matFF = matFF;
sysA.rhsI = rhsI;
sysA.rhsG = rhsG;
sysA.rhsF = rhsF;

% Full system
sysA.matA = [ matII matIG matIF; matGI matGG matGF; matFI matFG matFF ];
sysA.rhsA = [ rhsI ; rhsG ; rhsF];
spy(sysA.matA)

% % Reduced system
% sysA.matS = matGG - [matGI matGF]*([matII matIF ; matFI matFF]\[matIG ; matFG]);
% sysA.rhsS = rhsG - [matGI matGF]*([matII matIF ; matFI matFF]\[rhsI ; rhsF]);

% % Preconditionning
% if (PREC == 1)
%     sysA.matP = matGG;
%     sysA.matPinv = matGGinv;
% else
%     sysA.matP = 1;
%     sysA.matPinv = 1;
% end

% Compute solution
% solG = sysA.matS\sysA.rhsS;
% [solI ; solF] = [matII matIF ; matFI matFF]\([rhsI ; rhsF] - [matIG ; matFG]*solG);

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