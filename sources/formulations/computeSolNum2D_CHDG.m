function [solI, sysA, condLoc] = computeSolNum2D_CHDG(mesh, dofm, tau, BASIS, PREC)

global k edgTagToBC
global LdomX LdomY LpmlX LpmlY Rdom Rpml

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

    % Source terms
    rhsQ = mySourceVolume(xQ, yQ);

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

        % Mass matrices (physical space)
        matM_IIel = transpose(shapePhyQ) * (weightsQ .* shapePhyQ);
        matM_IGel = transpose(shapePhyQ) * (weightsQ .* shapeAuxQ);
        matM_GIel = transpose(shapeAuxQ) * (weightsQ .* shapePhyQ);
        matM_GGel = transpose(shapeAuxQ) * (weightsQ .* shapeAuxQ);

        % -----------------------------------------------------------------
        % Physical equations
        % -----------------------------------------------------------------

        % Local ID for interior unknowns and incoming characteristics
        idLocP = 0*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocU = 1*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocV = 2*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocG = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;

        % Element matrices (local element-wise system)
        matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 0.5*tau           * matM_IIel;
        matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + 0.5     * nx      * matM_IIel;
        matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + 0.5     * ny      * matM_IIel;
        matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 0.5               * matM_IGel;

        matIIel(idLocU,idLocP) = matIIel(idLocU,idLocP) + 0.5          * nx * matM_IIel;
        matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + 0.5/tau * nx * nx * matM_IIel;
        matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + 0.5/tau * nx * ny * matM_IIel;
        matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + 0.5/tau      * nx * matM_IGel;

        matIIel(idLocV,idLocP) = matIIel(idLocV,idLocP) + 0.5          * ny * matM_IIel;
        matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + 0.5/tau * nx * ny * matM_IIel;
        matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + 0.5/tau * ny * ny * matM_IIel;
        matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + 0.5/tau      * ny * matM_IGel;

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
            [souQ, ~, ~, souUxQ, souUyQ] = mySourceSurface(xQ, yQ);
            rhsPel = transpose(shapeAuxQ) * (weightsQ .* souQ);
            rhsNUel = transpose(shapeAuxQ) * (weightsQ .* (nx.*souUxQ + ny.*souUyQ));

            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));

            % Elemental matrices and RHS vectors (boundary conditions)
            matGGel = transpose(shapeAuxQ) * (weightsQ .* shapeAuxQ);
            matGIel = zeros(dofm.numDofPerLIN,3*dofm.numDofPerTRI);
            rhsGel = zeros(dofm.numDofPerLIN,1);
            switch BC
                case 'DIR0'
                    matGIel = [+tau*matM_GIel, +nx*matM_GIel, +ny*matM_GIel];
                case 'DIR'
                    matGIel = [+tau*matM_GIel, +nx*matM_GIel, +ny*matM_GIel];
                    rhsGel = +2*tau*rhsPel;
                case 'NEU0'
                    matGIel = [-tau*matM_GIel, -nx*matM_GIel, -ny*matM_GIel];
                case 'NEU'
                    matGIel = [-tau*matM_GIel, -nx*matM_GIel, -ny*matM_GIel];
                    rhsGel = -2*rhsNUel;
                case 'ABC'
                    matGIel = [+tau*matM_GIel, +nx*matM_GIel, +ny*matM_GIel] * (1-tau)/(1+tau);
                case 'ROB'
                    matGIel = [+tau*matM_GIel, +nx*matM_GIel, +ny*matM_GIel] * (1-tau)/(1+tau);
                    rhsGel = +(rhsPel - rhsNUel) * (2*tau)/(1+tau);
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
