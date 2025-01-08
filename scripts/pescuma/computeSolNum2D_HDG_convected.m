% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave, Simone Pescuma

function [solI, sysA] = computeSolNum2D_HDG_convected(mesh, dofm, PREC)

global edgTagToBC;
global omega c rho v0;

% -------------------------------------------------------------------------
% Build system
% -------------------------------------------------------------------------

numDofTRI = dofm.numDofTRI;
numDofFAC = dofm.numDofFAC;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);

% Shape functions and derivatives (reference space)
shapeLinRefQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeTriRefQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeTriRefDuQ, shapeTriRefDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% Global matrices
matIIx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIGx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIGy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIGv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN);
matIHx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN); %% new
matIHy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN); %% new
matIHv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerLIN); %% new
matGIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGGx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGGv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matGHx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN); %% new
matGHy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN); %% new
matGHv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN); %% new
matHIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN); %% new
matHIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN); %% new
matHIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN); %% new
matHGx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN); %% new
matHGy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN); %% new
matHGv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN); %% new
matHHx = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN); %% new
matHHy = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN); %% new
matHHv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN); %% new
matIIvInv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matGGvInv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matHHvInv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN); %% new

% Global RHS vectors
rhsI = zeros(3*numDofTRI,1);
rhsG = zeros(numDofFAC,1);
rhsH = zeros(numDofFAC,1); %% new

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
    shapeTriQ = shapeTriRefQ * orientation;
    shapeDxQ = (shapeTriRefDuQ * Jdudx(1,1) + shapeTriRefDvQ * Jdudx(2,1)) * orientation;
    shapeDyQ = (shapeTriRefDuQ * Jdudx(1,2) + shapeTriRefDvQ * Jdudx(2,2)) * orientation;

    % Source terms
    rhsQ = mySourceVolume(xQ, yQ);

    % Elemental matrices/vectors
    matMel = shapeTriQ' * (weightsTriQ .* shapeTriQ) * detJdxdu;
    matDXel = shapeDxQ' * (weightsTriQ .* shapeTriQ) * detJdxdu;
    matDYel = shapeDyQ' * (weightsTriQ .* shapeTriQ) * detJdxdu;
    rhsPel = shapeTriQ' * (weightsTriQ .* rhsQ) * detJdxdu;

    matIIel = [
        -1i*omega/(rho*c^2)*matMel-1/(rho*c^2)*(matDXel*v0(1)+matDYel*v0(2))   -matDXel                                                 -matDYel                         ;
        -matDXel                                                               -1i*omega*rho*matMel-rho*(matDXel*v0(1)+matDYel*v0(2))   zeros(numDofPerTRI,numDofPerTRI) ;
        -matDYel                                                               zeros(numDofPerTRI,numDofPerTRI)                         -1i*omega*rho*matMel-rho*(matDXel*v0(1)+matDYel*v0(2))];

    rhsIel = [
        -1/(1i*omega*rho)*rhsPel ;
        zeros(numDofPerTRI,1) ;
        zeros(numDofPerTRI,1) ];

    % ---------------------------------------------------------------------
    % Surface terms
    % ---------------------------------------------------------------------

    matIGel = zeros(3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);
    matIHel = zeros(3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);  %% new

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
        shapeLinQ = shapeLinRefQ * orientation;

        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);

        % Exterior tangent
        tx = -ny;
        ty = nx;

        v0n = v0(1)*nx+v0(2)*ny;

        % Mass matrices (physical space)
        matMel = shapeLinQ' * (weightsLinQ .* shapeLinQ) * Jdxdu;

        % -----------------------------------------------------------------
        % Physical equations
        % -----------------------------------------------------------------

        % Local ID for interior and auxiliary unknowns
        idLocP = 0*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocU = 1*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocV = 2*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocG = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;
        idLocH = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN; %% new

        % Element matrices (local element-wise system)
        matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 1/(rho*c)                      * matMel;  %%
        matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) +                             nx * matMel;  %%
        matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) +                             ny * matMel;  %%

        matIIel(idLocU,idLocP) = matIIel(idLocU,idLocP) + v0n / c                   * nx * matMel;  %% 
        matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + v0n*rho              * nx * nx * matMel;  %%
        matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + v0n*rho              * ny * nx * matMel;  %%

        matIIel(idLocV,idLocP) = matIIel(idLocV,idLocP) + v0n / c                   * ny * matMel;  %%
        matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + v0n*rho              * nx * ny * matMel;  %%
        matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + v0n*rho              * ny * ny * matMel;  %%

        matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) + (v0n-c) / (rho*c^2)            * matMel;  %%
        matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + (c-v0n) / c               * nx * matMel;  %%
        matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + (c-v0n) / c               * ny * matMel;  %%

        matIHel(idLocP,idLocH) = matIHel(idLocP,idLocH) + 0*matMel;  %%
        matIHel(idLocU,idLocH) = matIHel(idLocU,idLocH) + v0n*rho * tx * matMel;  %%
        matIHel(idLocV,idLocH) = matIHel(idLocV,idLocH) + v0n*rho * ty * matMel;  %%

        % -----------------------------------------------------------------
        % Auxiliary equations
        % -----------------------------------------------------------------

        matGHel = 0*matMel;  %%%%
        matGIel = zeros(dofm.numDofPerLIN,3*dofm.numDofPerLIN);  %%%%
        rhsGel  = zeros(dofm.numDofPerLIN,1);  %%%%
        matHHel = matMel;  %%%%
        matHGel = 0*matMel;  %%%%%
        matHIel = zeros(dofm.numDofPerLIN,3*dofm.numDofPerLIN);  %%%%
        rhsHel  = zeros(dofm.numDofPerLIN,1);  %%%%

        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);

        if (triNeigh > 0)

            matGGel = 0.5 * matMel;
            matGIel = - 0.5 * [matMel, rho*c*nx*matMel, rho*c*ny*matMel];

            % Global ID for auxiliary and interior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            idGloH = dofm.locToGloFAC(tri,idLocH);
            idGloP = 0*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloU = 1*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloV = 2*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloI = [idGloP idGloU idGloV];

            % Assembling
            matGIx(idGloG,:) = idGloG'*ones(1,size(idGloI,2));
            matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
            matGHx(idGloG,:) = idGloG'*ones(1,size(idGloH,2)); %% new
            matGIy(idGloG,:) = ones(size(idGloG,2),1)*idGloI;
            matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
            matGHy(idGloG,:) = ones(size(idGloG,2),1)*idGloH;  %% new
            matGIv(idGloG,:) = matGIel;
            matGGv(idGloG,:) = matGGel;
            matGHv(idGloG,:) = matGHel;  %% new
            matGGvInv(idGloG,:) = inv(matGGel);
            rhsG(idGloG) = rhsGel;

            matHIel = - [0*matMel, tx*matMel, ty*matMel];

            if v0n<0
                % Global ID for auxiliary and exterior unknowns
                idGloG = dofm.locToGloFAC(tri,idLocG);
                idGloH = dofm.locToGloFAC(tri,idLocH);
                dofExt = dofm.locFacNeigh(facNeigh,:);
                idGloP = 0*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
                idGloU = 1*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
                idGloV = 2*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
                idGloI = [idGloP idGloU idGloV];
            end

            matHIx(idGloH,:) = idGloH'*ones(1,size(idGloI,2));  %% new
            matHGx(idGloH,:) = idGloH'*ones(1,size(idGloG,2));  %% new
            matHHx(idGloH,:) = idGloH'*ones(1,size(idGloH,2));  %% new
            matHIy(idGloH,:) = ones(size(idGloH,2),1)*idGloI;  %% new
            matHGy(idGloH,:) = ones(size(idGloH,2),1)*idGloG;  %% new
            matHHy(idGloH,:) = ones(size(idGloH,2),1)*idGloH;  %% new
            matHIv(idGloH,:) = matHIel;  %% new
            matHGv(idGloH,:) = matHGel;  %% new
            matHHv(idGloH,:) = matHHel;  %% new
            matHHvInv(idGloH,:) = inv(matHHel); %% new
            rhsH(idGloH) = rhsHel;  %% new

        else

            matGGel = matMel;

            % Source terms
            [solQ, ~, ~, solVxQ, solVyQ] = mySourceSurface(xQ,yQ);
            rhsPel = shapeLinQ' * (weightsLinQ .* solQ) * Jdxdu;
            rhsUel = shapeLinQ' * (weightsLinQ .* solVxQ) * Jdxdu;
            rhsVel = shapeLinQ' * (weightsLinQ .* solVyQ) * Jdxdu;

            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));

            if v0n>0
                matHIel = - [0*matMel, tx*matMel, ty*matMel];
            else
                rhsHel = tx * rhsUel + ty * rhsVel;
            end

            % Elemental matrices/vectors (boundary condition)
            switch BC
                case 'DIR'
                    rhsGel = rhsPel; 
                case 'NEU'
                    matGIel = - [matMel, rho*c*nx*matMel, rho*c*ny*matMel];
                    rhsGel = -rho*c*(nx*rhsUel + ny*rhsVel);
                case 'NEU0'
                    matGIel = - [matMel, rho*c*nx*matMel, rho*c*ny*matMel];
                    rhsGel = 0 * rhsUel + 0 * rhsVel;
                case 'ABC'
                    matGIel = -1/2 * [matMel, rho*c*nx*matMel, rho*c*ny*matMel];
                    rhsGel = 0 * rhsPel;
                case 'ROB'
                    matGIel = -1/2 * [matMel, rho*c*nx*matMel, rho*c*ny*matMel];
                    rhsGel = (rhsPel - rho*c * (nx*rhsUel + ny*rhsVel)) / 2;
                otherwise
                    error('BAD BOUNDARY CONDITION.');
            end

            % Global ID for auxiliary unknowns and interior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            idGloH = dofm.locToGloFAC(tri,idLocH);  %% new
            idGloP = 0*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloU = 1*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloV = 2*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloI = [idGloP idGloU idGloV];

            % Assembling
            matGIx(idGloG,:) = idGloG'*ones(1,size(idGloI,2));
            matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
            matGHx(idGloG,:) = idGloG'*ones(1,size(idGloH,2)); %% new
            matGIy(idGloG,:) = ones(size(idGloG,2),1)*idGloI;
            matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
            matGHy(idGloG,:) = ones(size(idGloG,2),1)*idGloH;  %% new
            matGIv(idGloG,:) = matGIel;
            matGGv(idGloG,:) = matGGel;
            matGHv(idGloG,:) = matGHel;  %% new
            matGGvInv(idGloG,:) = inv(matGGel);
            rhsG(idGloG) = rhsGel;
            
            matHIx(idGloH,:) = idGloH'*ones(1,size(idGloI,2));  %% new
            matHGx(idGloH,:) = idGloH'*ones(1,size(idGloG,2));  %% new
            matHHx(idGloH,:) = idGloH'*ones(1,size(idGloH,2));  %% new
            matHIy(idGloH,:) = ones(size(idGloH,2),1)*idGloI;  %% new
            matHGy(idGloH,:) = ones(size(idGloH,2),1)*idGloG;  %% new
            matHHy(idGloH,:) = ones(size(idGloH,2),1)*idGloH;  %% new
            matHIv(idGloH,:) = matHIel;  %% new
            matHGv(idGloH,:) = matHGel;  %% new
            matHHv(idGloH,:) = matHHel;  %% new
            matHHvInv(idGloH,:) = inv(matHHel); %% new
            rhsH(idGloH) = rhsHel;  %% new

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
    dofGloH = dofm.locToGloFAC(tri,:);  %% new

    % Assembling
    idTRI = (tri-1)*3*dofm.numDofPerTRI + (1:3*dofm.numDofPerTRI);
    matIIx(idTRI,:) = dofGloI'*ones(1,size(dofGloI,2));
    matIGx(idTRI,:) = dofGloI'*ones(1,size(dofGloG,2));
    matIHx(idTRI,:) = dofGloI'*ones(1,size(dofGloH,2));  %% new
    matIIy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloI;
    matIGy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloG;
    matIHy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloH;  %% new
    matIIv(idTRI,:) = matIIel;
    matIGv(idTRI,:) = matIGel;
    matIHv(idTRI,:) = matIHel;  %% new
    matIIvInv(idTRI,:) = inv(matIIel);
    rhsI(dofGloI) = rhsIel;

end

% Sparse memory storage
matII = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofFAC);
matIH = sparse(matIHx, matIHy, matIHv, 3*numDofTRI, numDofFAC);  %% new
matGI = sparse(matGIx, matGIy, matGIv, numDofFAC, 3*numDofTRI);
matGG = sparse(matGGx, matGGy, matGGv, numDofFAC, numDofFAC);
matGH = sparse(numDofFAC, numDofFAC);  %% new
matHI = sparse(matHIx, matHIy, matHIv, numDofFAC, 3*numDofTRI);  %% new
matHG = sparse(numDofFAC, numDofFAC);  %% new
matHH = sparse(matHHx, matHHy, matHHv, numDofFAC, numDofFAC);  %% new
matIIinv = sparse(matIIx, matIIy, matIIvInv, 3*numDofTRI, 3*numDofTRI);
matGGinv = sparse(matGGx, matGGy, matGGvInv, numDofFAC, numDofFAC);
matHHinv = sparse(matHHx, matHHy, matHHvInv, numDofFAC, numDofFAC);

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

% Matrix partition
sysA.matII = matII;
sysA.matIG = [matIG matIH];
sysA.matGI = [matGI; matHI];
sysA.matGG = [matGG matGH; matHG matHH];
sysA.rhsI = rhsI;
sysA.rhsG = [rhsG; rhsH];
sysA.matIIinv = matIIinv;
sysA.matGGinv = blkdiag(matGGinv, matHHinv);

% Full system
sysA.matA = [ matII matIG matIH; matGI matGG matGH; matHI matHG matHH];
sysA.rhsA = [ rhsI ; rhsG ; rhsH ];

% Preconditionning matrix
if(PREC == 1)
    sysA.matP = blkdiag(matGG, matHH);
    sysA.matPinv = blkdiag(matGGinv, matHHinv);
else
    sysA.matP = eye(size(matGG,1)+size(matHH,1));
    sysA.matPinv = eye(size(matGG,1)+size(matHH,1));
end

% Reduced system
sysA.matS = [matGG-matGI*(matIIinv*matIG) -matGI*(matIIinv*matIH);
             -matHI*(matIIinv*matIG)      matHH-matHI*(matIIinv*matIH)];
sysA.rhsS = [rhsG-matGI*(matIIinv*rhsI); rhsH-matHI*(matIIinv*rhsI)];

% Physical system
sysA.matPhy = matII - matIG*(matGGinv*matGI) - matIH*(matHHinv*matHI);
sysA.rhsPhy = rhsI - matIG*(matGGinv*rhsG) - matIH*(matHHinv*rhsH);

% Compute solution
solX = sysA.matS\sysA.rhsS;
L_G = size(rhsG-matGI*(matIIinv*rhsI),1);
L_H = size(rhsH-matHI*(matIIinv*rhsI),1);
solG = solX(1:L_G);
solH = solX(1+L_G:L_G+L_H);
solI = matIIinv*(rhsI-matIG*solG-matIH*solH);

% A = [matII matIG;
%      matGI  matGG];
% rhs = [rhsI;
%        rhsG];
% 
% solX = A\rhs;

end