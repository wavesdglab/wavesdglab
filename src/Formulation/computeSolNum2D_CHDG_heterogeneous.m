% Contributors: Simone Pescuma, Axel Modave

function [solI, sysA] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, FLUX, PREC)

% -------------------------------------------------------------------------
% Build system
% -------------------------------------------------------------------------

disp('--- Build system ---');

global eta k omega edgTagToBC

numDofTRI = dofm.numDofTRI;
numDofFAC = dofm.numDofFAC;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);

% Shape functions and derivatives (reference space)
shapeLinRefQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeLinRefDuQ = functionsShapeDerLIN(uLinQ, dofm.degree);
shapeTriRefQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeTriRefDuQ, shapeTriRefDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

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
matPPv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
matIIvInv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matGGvInv = zeros(mesh.numTri*3*dofm.numDofPerLIN, dofm.numDofPerLIN);
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
        -1i*k(tri)/eta(tri)*matMel  -matDXel                          -matDYel                         ;
        -matDXel                    -1i*k(tri)*eta(tri)*matMel        zeros(numDofPerTRI,numDofPerTRI) ;
        -matDYel                    zeros(numDofPerTRI,numDofPerTRI)  -1i*k(tri)*eta(tri)*matMel       ];

    rhsIel = [
        -1/(1i*k(tri)*eta(tri))*rhsPel ;
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
        if(n1(fac) > n2(fac))
            orientation(3:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        orientation = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation);

        % Shape functions (physical space)
        shapeLinQ = shapeLinRefQ * orientation;
        shapeLinDsQ = shapeLinRefDuQ * orientation;

        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);

        % Mass matrices (physical space)
        matMel = shapeLinQ' * (weightsLinQ .* shapeLinQ) * Jdxdu;
        matKel = shapeLinDsQ' * (weightsLinQ .* shapeLinDsQ) / Jdxdu;

        % -----------------------------------------------------------------
        % Physical equations
        % -----------------------------------------------------------------

        triNeigh = mesh.mapTriToTri(tri,fac);
        if (triNeigh>0)
            etaNeigh = eta(triNeigh);
            kNeigh = k(triNeigh);
        else
            etaNeigh = eta(tri);
            kNeigh = k(tri);
        end

        % Local ID for interior and auxiliary unknowns
        idLocP = 0*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocU = 1*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocV = 2*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocG = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;

        % Element matrices (local element-wise system)
        switch FLUX

            case 'UPW'

                matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 1/(eta(tri) + etaNeigh)                      * matMel;
                matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + eta(tri)/(eta(tri) + etaNeigh)               * nx * matMel;
                matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + eta(tri)/(eta(tri) + etaNeigh)               * ny * matMel;
                matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 1/(eta(tri) + etaNeigh)                      * matMel;

                matIIel(idLocU,idLocP) = matIIel(idLocU,idLocP) + etaNeigh/(eta(tri) + etaNeigh)          * nx * matMel;
                matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + etaNeigh*eta(tri)/(eta(tri) + etaNeigh) * nx * nx * matMel;
                matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + etaNeigh*eta(tri)/(eta(tri) + etaNeigh) * nx * ny * matMel;
                matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + eta(tri)/(eta(tri) + etaNeigh)               * nx * matMel;

                matIIel(idLocV,idLocP) = matIIel(idLocV,idLocP) + etaNeigh/(eta(tri) + etaNeigh)          * ny * matMel;
                matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + etaNeigh*eta(tri)/(eta(tri) + etaNeigh) * nx * ny * matMel;
                matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + etaNeigh*eta(tri)/(eta(tri) + etaNeigh) * ny * ny * matMel;
                matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + eta(tri)/(eta(tri) + etaNeigh)               * ny * matMel;

            case 'UPW2'

                matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 0.5/(eta(tri))         * matMel;
                matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + 0.5               * nx * matMel;
                matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + 0.5               * ny * matMel;
                matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 0.5/etaNeigh           * matMel;

                matIIel(idLocU,idLocP) = matIIel(idLocU,idLocP) + 0.5               * nx * matMel;
                matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + 0.5*eta(tri) * nx * nx * matMel;
                matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + 0.5*eta(tri) * ny * nx * matMel;
                matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + 0.5               * nx * matMel;

                matIIel(idLocV,idLocP) = matIIel(idLocV,idLocP) + 0.5               * ny * matMel;
                matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + 0.5*eta(tri) * nx * ny * matMel;
                matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + 0.5*eta(tri) * ny * ny * matMel;
                matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + 0.5               * ny * matMel;

            case 'SYM'

                etaF = sqrt(eta(tri)*etaNeigh);

                matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 0.5/etaF          * matMel;
                matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + 0.5     * nx      * matMel;
                matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + 0.5     * ny      * matMel;
                matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 0.5               * matMel;

                matIIel(idLocU,idLocP) = matIIel(idLocU,idLocP) + 0.5          * nx * matMel;
                matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + 0.5*etaF* nx * nx * matMel;
                matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + 0.5*etaF* nx * ny * matMel;
                matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + 0.5*etaF     * nx * matMel;

                matIIel(idLocV,idLocP) = matIIel(idLocV,idLocP) + 0.5          * ny * matMel;
                matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + 0.5*etaF* nx * ny * matMel;
                matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + 0.5*etaF* ny * ny * matMel;
                matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + 0.5*etaF     * ny * matMel;

            case 'SYM2'

                etaF = sqrt(eta(tri)*etaNeigh);
                kF = sqrt(k(tri)*kNeigh);
                matEtaF = etaF * matMel / (matMel + 0.5/(kF^2) * matKel);

                matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 0.5 *         (matEtaF\matMel);
                matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + 0.5         * nx      * matMel;
                matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + 0.5         * ny      * matMel;
                matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 0.5                   * matMel;

                matIIel(idLocU,idLocP) = matIIel(idLocU,idLocP) + 0.5              * nx * matMel;
                matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + 0.5*matEtaF * nx * nx * matMel;
                matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + 0.5*matEtaF * nx * ny * matMel;
                matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + 0.5*matEtaF      * nx * matMel;

                matIIel(idLocV,idLocP) = matIIel(idLocV,idLocP) + 0.5              * ny * matMel;
                matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + 0.5*matEtaF * nx * ny * matMel;
                matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + 0.5*matEtaF * ny * ny * matMel;
                matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + 0.5*matEtaF      * ny * matMel;

            otherwise
                error('BAD FLUX.');
        end

        % -----------------------------------------------------------------
        % Auxiliary equations
        % -----------------------------------------------------------------

        % Elemental matrices/vectors
        matGGel = matMel;
        matGIel = zeros(dofm.numDofPerLIN,3*dofm.numDofPerLIN);
        rhsGel = zeros(dofm.numDofPerLIN,1);

        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);

        if (triNeigh > 0)

            % Elemental matrices (interface condition)
            switch FLUX
                case {'UPW','UPW2'}
                    matGIel = [-matMel, eta(triNeigh)*nx*matMel, eta(triNeigh)*ny*matMel];
                case 'SYM'
                    matGIel = [-1/etaF*matMel, nx*matMel, ny*matMel];
                case 'SYM2'
                    matGIel = [-(matEtaF\matMel), nx*matMel, ny*matMel];
                otherwise
                    error('BAD FLUX.');
            end

            % Global ID for auxiliary and exterior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            dofExt = dofm.locFacNeigh(facNeigh,:);
            idGloP = 0*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idGloU = 1*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idGloV = 2*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idGloI = [idGloP idGloU idGloV];

        else

            % Source terms
            [solP, ~, ~, solU, solV] = mySourceSurface(xQ, yQ);
            rhsPel = shapeLinQ' * (weightsLinQ .* solP) * Jdxdu;
            rhsNUel = shapeLinQ' * (weightsLinQ .* (nx.*solU + ny.*solV)) * Jdxdu;

            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));

            % Elemental matrices/vectors (boundary condition)
            switch FLUX
                case {'UPW','UPW2'}
                    switch BC
                        case 'DIR0'
                            matGIel = [matMel, eta(tri)*nx*matMel, eta(tri)*ny*matMel];
                        case 'DIR'
                            matGIel = [matMel, eta(tri)*nx*matMel, eta(tri)*ny*matMel];
                            rhsGel = 2*rhsPel;
                        case 'NEU0'
                            matGIel = [-matMel, -eta(tri)*nx*matMel, -eta(tri)*ny*matMel];
                        case 'NEU'
                            matGIel = [-matMel, -eta(tri)*nx*matMel, -eta(tri)*ny*matMel];
                            rhsGel = -2*eta(tri)*rhsNUel;
                        case 'ABC'
                        case 'ROB'
                            rhsGel = rhsPel - etaNeigh*rhsNUel;
                        otherwise
                            error('BAD BOUNDARY CONDITION.');
                    end
                case 'SYM'
                    switch BC
                        case 'DIR0'
                            matGIel = [eta(tri)\matMel, nx*matMel, ny*matMel];
                        case 'DIR'
                            matGIel = [eta(tri)\matMel, nx*matMel, ny*matMel];
                            rhsGel = 2/eta(tri)*rhsPel;
                        case 'NEU0'
                            matGIel = [-eta(tri)\matMel, -nx*matMel, -ny*matMel];
                        case 'NEU'
                            matGIel = [-eta(tri)\matMel, -nx*matMel, -ny*matMel];
                            rhsGel = -2*rhsNUel;
                        case 'ABC'
                        case 'ROB'
                            rhsGel = (rhsPel - eta(tri)*rhsNUel)/eta(tri);
                        otherwise
                            error('BAD BOUNDARY CONDITION.');
                    end
                case 'SYM2'
                    switch BC
                        case 'DIR0'
                            matGIel = [matEtaF\matMel, nx*matMel, ny*matMel];
                        case 'DIR'
                            matGIel = [matEtaF\matMel, nx*matMel, ny*matMel];
                            rhsGel = 2*matEtaF\rhsPel;
                        case 'NEU0'
                            matGIel = [-matEtaF\matMel, -nx*matMel, -ny*matMel];
                        case 'NEU'
                            matGIel = [-matEtaF\matMel, -nx*matMel, -ny*matMel];
                            rhsGel = -2*rhsNUel;
                        case 'ABC'
                            matTmpP = sparse(eye(size(matEtaF,1))) + matEtaF/eta(tri);
                            matTmpM = sparse(eye(size(matEtaF,1))) - matEtaF/eta(tri);
                            matGIel = (matTmpM/matTmpP) * [-matEtaF\matMel, -nx*matMel, -ny*matMel];
                        case 'ROB'
                            matTmpP = sparse(eye(size(matEtaF,1))) + matEtaF/eta(tri);
                            matTmpM = sparse(eye(size(matEtaF,1))) - matEtaF/eta(tri);
                            matGIel = (matTmpM/matTmpP) * [-matEtaF\matMel, -nx*matMel, -ny*matMel];
                            rhsGel  = (2/eta(tri)) * inv(matTmpP) * (rhsPel - eta(tri)*rhsNUel);
                        otherwise
                            error('BAD BOUNDARY CONDITION.');
                    end
                otherwise
                    error('BAD FLUX.');
            end

            % Global ID for auxiliary unknowns and interior unknowns
            idGloG = dofm.locToGloFAC(tri,idLocG);
            idGloP = 0*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloU = 1*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloV = 2*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloI = [idGloP idGloU idGloV];
        end

        % Preconditionning matrix
        if(PREC == 1)
            switch FLUX
                case 'UPW'
                    matPPel = (eta(tri) + etaNeigh)*matMel;
                case 'UPW2'
                    matPPel = eta(tri)\matMel;
                case {'SYM','SYM2'}
                    matPPel = matMel;
                otherwise
                    error('BAD FLUX.');
            end
        else
            matPPel = matMel;
        end
        
        % Assembling
        matGIx(idGloG,:) = idGloG'*ones(1,size(idGloI,2));
        matGGx(idGloG,:) = idGloG'*ones(1,size(idGloG,2));
        matGIy(idGloG,:) = ones(size(idGloG,2),1)*idGloI;
        matGGy(idGloG,:) = ones(size(idGloG,2),1)*idGloG;
        matGIv(idGloG,:) = matGIel;
        matGGv(idGloG,:) = matGGel;
        matPPv(idGloG,:) = matPPel;
        matGGvInv(idGloG,:) = inv(matGGel);
        matPPvInv(idGloG,:) = inv(matPPel);
        rhsG(idGloG) = rhsGel;
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
% Point source
% -------------------------------------------------------------------------

global pntSouTag pntSouVal
if(~isempty(pntSouTag))
    TOT = 0;
    SouEl = 0;
    for tri=1:mesh.numTri
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
    for ind=1:TOT
        dofSou = dofm.numDofPerTRI * (SouEl(1,ind)-1) + SouEl(2,ind);
        rhsI(dofSou) = rhsI(dofSou) - 1/(1i*omega) * pntSouVal / TOT;
    end
end

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

disp('--- Solve system ---');

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
sysA.matP = matPP;
sysA.matPinv = matPPinv;
%sysA.matP = matGG;
%sysA.matPinv = matGGinv;

% Compute solution
solG = sysA.matS\sysA.rhsS;
solI = matIIinv*(rhsI-matIG*solG);

end
