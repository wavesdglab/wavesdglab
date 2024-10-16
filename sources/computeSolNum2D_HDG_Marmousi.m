% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solI, sysA] = computeSolNum2D_HDG_Marmousi(mesh, dofm, PREC)

global omega edgTagToBC
global etaArray kArray
global pntSouTag pntSouVal

eta = etaArray;
k = kArray;

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
[shapeTriDuQ, shapeTriDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% Global matrices
matIIx = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);
matIIy = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);
matIIv = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);
matIGx = zeros(3*dofm.numDofPerTRI,mesh.numTri*3*dofm.numDofPerLIN);
matIGy = zeros(3*dofm.numDofPerTRI,mesh.numTri*3*dofm.numDofPerLIN);
matIGv = zeros(3*dofm.numDofPerTRI,mesh.numTri*3*dofm.numDofPerLIN);
matGIx = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);
matGIy = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);
matGIv = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);
matGGx = zeros(numDofLIN,dofm.numDofPerLIN);
matGGy = zeros(numDofLIN,dofm.numDofPerLIN);
matGGv = zeros(numDofLIN,dofm.numDofPerLIN);
matIIvInv = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);
matGGvInv = zeros(numDofLIN,dofm.numDofPerLIN);

% Global RHS vectors
rhsI = zeros(3*numDofTRI,1);
rhsG = zeros(numDofLIN,1);

TOT = 0;
SouEl = 0;

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

    if(~isempty(pntSouTag))
        vertSou = mesh.mapPntToVer(mesh.tagPntFile == pntSouTag);
        for pos = 1:3
            if(mesh.mapTriToVer(tri,pos) == vertSou)
                TOT = TOT+1;
                SouEl(1,TOT) = tri;
                SouEl(2,TOT) = pos;
                SouEl(3,TOT) = 0;
                if (V1(1,2) == 0 && V2(1,2) == 0) || (V1(1,2) == 0 && V3(1,2) == 0) || (V2(1,2) == 0 && V3(1,2) == 0)
                    SouEl(3,TOT) = 1;
                end
            end
        end
    end

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
    rhsQ = mySourceVolume(xQ, yQ);
    
    % Elemental matrices and RHS vectors
    matMel = shapePhyQ' * (weightsTriQ .* shapePhyQ) * detJdxdu;
    matDXel = shapeDxQ' * (weightsTriQ .* shapePhyQ) * detJdxdu;
    matDYel = shapeDyQ' * (weightsTriQ .* shapePhyQ) * detJdxdu;
    rhsPel = shapePhyQ' * (weightsTriQ .* rhsQ) * detJdxdu;
    
    matIIel = [
        -1i*k(tri)/eta(tri)*matMel  -matDXel                          -matDYel                         ;
        -matDXel                    -1i*k(tri)*eta(tri)*matMel        zeros(numDofPerTRI,numDofPerTRI) ;
        -matDYel                    zeros(numDofPerTRI,numDofPerTRI)  -1i*k(tri)*eta(tri)*matMel       ];
   
    rhsIel = [
        -1/(1i*k(tri)*eta(tri))*rhsPel  ;
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
        
        triNeigh = mesh.mapTriToTri(tri,fac);

        if (triNeigh>0)
            etaF = sqrt(eta(tri)*eta(triNeigh));
            kF = sqrt(k(tri)*k(triNeigh));
        else
            etaF = eta(tri);
            kF = k(tri);
        end

        % Mapping
        V1 = mesh.coord(n1(fac),:);
        V2 = mesh.coord(n2(fac),:);
        [xQ, yQ] = locToGloLIN(uLinQ, V1, V2);
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
        
        % Mass matrices (physical space)
        matM_IIel = shapePhyQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_IGel = shapePhyQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        matM_GIel = shapeAuxQ' * (weightsLinQ .* shapePhyQ) * Jdxdu;
        matM_GGel = shapeAuxQ' * (weightsLinQ .* shapeAuxQ) * Jdxdu;
        
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
        matIIel(dofLocP,dofLocP) = matIIel(dofLocP,dofLocP) + 1 / etaF * matM_IIel;
        matIIel(dofLocP,dofLocU) = matIIel(dofLocP,dofLocU) +      nx * matM_IIel;
        matIIel(dofLocP,dofLocV) = matIIel(dofLocP,dofLocV) +      ny * matM_IIel;
        matIGel = zeros(3*dofm.numDofPerTRI,dofm.numDofPerLIN);
        matIGel(dofLocP,:) = - 1 / etaF * matM_IGel;
        matIGel(dofLocU,:) =       nx  * matM_IGel;
        matIGel(dofLocV,:) =       ny  * matM_IGel;
        
        % -----------------------------------------------------------------
        % Auxiliary equations
        % -----------------------------------------------------------------
        
        matGGel = zeros(dofm.numDofPerLIN,dofm.numDofPerLIN);
        matGIel = zeros(dofm.numDofPerLIN,3*dofm.numDofPerLIN);
        rhsGel  = zeros(dofm.numDofPerLIN,1);
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        
        if (triNeigh > 0)
            
            % Elemental matrices (interface condition)
            matGGel = 0.5 * matM_GGel;
            matGIel = - etaF / 2 * [1/etaF * matM_GIel, nx*matM_GIel, ny*matM_GIel];
            
        else
            
            % Source terms
            [solQ, solDxQ, solDyQ, ~, ~] = mySol(xQ, yQ);
            rhsPel = shapeAuxQ' * (weightsLinQ .* solQ) * Jdxdu;
            rhsUel = shapeAuxQ' * (weightsLinQ .* solDxQ) * Jdxdu / (1i*kF*etaF);
            rhsVel = shapeAuxQ' * (weightsLinQ .* solDyQ) * Jdxdu / (1i*kF*etaF);
            
            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));
            
            % Elemental matrices and RHS vectors (boundary conditions)
            matGGel = matM_GGel;
            switch BC
                case 'DIR0'
                case 'DIR'
                    rhsGel = rhsPel;
                case 'NEU'
                    matGIel = - [matM_GIel, etaF*nx*matM_GIel, etaF*ny*matM_GIel];
                    rhsGel = -etaF*(nx*rhsUel + ny*rhsVel);
                case 'ABC'
                    matGIel = -1/2 * [matM_GIel, etaF*nx*matM_GIel, etaF*ny*matM_GIel];
%                     rhsGel = (rhsPel - etaF * (nx*rhsUel + ny*rhsVel)) / 2;
                    rhsGel = 0*(rhsPel - etaF * (nx*rhsUel + ny*rhsVel)) / 2;
                case 'ROB'
                    matGIel = -1/2 * [matM_GIel, etaF*nx*matM_GIel, etaF*ny*matM_GIel];
                    rhsGel = (rhsPel - etaF * (nx*rhsUel + ny*rhsVel)) / 2;    
                otherwise
                    error('BAD BOUNDARY CONDITION.');
            end
        end
        
        % Global ID for edge unknowns
        edgGlo = abs(mesh.mapTriToEdg(tri,fac));
        dofGloG = dofm.locToGloLIN(edgGlo,:);
        dofLocG = 1:dofm.numDofPerLIN;
        if(mesh.mapTriToEdg(tri,fac) < 0)
            dofLocG(1) = 2;
            dofLocG(2) = 1;
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

        if triNeigh > 0
            matGGvInv(dofGloG,:) = 1/2*inv(matGGel(dofLocG,dofLocG));  % NEW
        else
            matGGvInv(dofGloG,:) = matGGvInv(dofGloG,:) + inv(matGGel(dofLocG,dofLocG));  % NEW
        end
        % it works, but why?

%         matGGvInv(dofGloG,:) = matGGvInv(dofGloG,:) + inv(matGGel(dofLocG,dofLocG));  % NEW
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
    
end

matII    = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG    = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofLIN);
matGI    = sparse(matGIx, matGIy, matGIv, numDofLIN, 3*numDofTRI);
matGG    = sparse(matGGx, matGGy, matGGv, numDofLIN, numDofLIN);
matIIinv = sparse(matIIx, matIIy, matIIvInv, 3*numDofTRI, 3*numDofTRI);
matGGinv = sparse(matGGx, matGGy, matGGvInv, numDofLIN, numDofLIN);

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

% THE GOOD ONE!
if(~isempty(pntSouTag))
    for ind=1:TOT
%         if (SouEl(3,ind) == 1)
            dofSou = dofm.numDofPerTRI * (SouEl(1,ind)-1) + SouEl(2,ind);
            % SouEl(1,ind) = tri; SouEl(2,ind) = DOF associated to the node for pressure
            rhsI(dofSou) = rhsI(dofSou) - 1/(1i*omega) * pntSouVal / TOT;
%         end
    end
end

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
disp('--- Full system ---');
tic
sysA.matA = [ matII matIG ; matGI matGG ];
sysA.rhsA = [ rhsI ; rhsG ];
toc

% Reduced system
disp('--- Reduced system ---');
tic
sysA.matS = matGG - matGI*(matIIinv*matIG);
sysA.rhsS = rhsG - matGI*(matIIinv*rhsI);
toc

% Physical system
disp('--- Physical system ---');
tic
sysA.matPhy = matII - matIG*(sysA.matGGinv*matGI);
sysA.rhsPhy = rhsI - matIG*(sysA.matGGinv*rhsG);
toc

% Preconditioning
disp('--- Preconditioning ---');
tic
if (PREC == 1)
    sysA.matP = matGG;
    sysA.matPinv = sysA.matGGinv;
else
    sysA.matP = 1;
    sysA.matPinv = 1;
end
toc

% Compute solution
disp('--- Compute direct solution ---');
tic
solG = sysA.matS\sysA.rhsS;
solI = matIIinv*(rhsI-matIG*solG);
toc

end