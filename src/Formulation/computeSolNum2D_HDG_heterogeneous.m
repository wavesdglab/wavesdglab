% Contributors: Simone Pescuma, Axel Modave

function [solI, sysA] = computeSolNum2D_HDG_heterogeneous(mesh, dofm, FLUX)

% -------------------------------------------------------------------------
% Build system
% -------------------------------------------------------------------------

disp('--- Build system ---');

global eta k omega edgTagToBC

numDofTRI = dofm.numDofTRI;
numDofLIN = dofm.numDofLIN;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);

% Shape functions and derivatives (reference space)
shapeRefLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeRefTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeTriDuQ, shapeTriDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% Global matrices
matIIx = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIy = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIIv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matIGx = zeros(3*dofm.numDofPerTRI, mesh.numTri*3*dofm.numDofPerLIN);
matIGy = zeros(3*dofm.numDofPerTRI, mesh.numTri*3*dofm.numDofPerLIN);
matIGv = zeros(3*dofm.numDofPerTRI, mesh.numTri*3*dofm.numDofPerLIN);
matGIx = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIy = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGIv = zeros(mesh.numTri*3*dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
matGGx = zeros(numDofLIN, dofm.numDofPerLIN);
matGGy = zeros(numDofLIN, dofm.numDofPerLIN);
matGGv = zeros(numDofLIN, dofm.numDofPerLIN);
matIIvInv = zeros(mesh.numTri*3*dofm.numDofPerTRI, 3*dofm.numDofPerTRI);
matGGvInv = zeros(numDofLIN, dofm.numDofPerLIN);

% Global RHS vectors
rhsI = zeros(3*numDofTRI,1);
rhsG = zeros(numDofLIN,1);

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
    shapePhyQ = shapeRefTriQ * orientation;
    shapeDxQ = (shapeTriDuQ * Jdudx(1,1) + shapeTriDvQ * Jdudx(2,1)) * orientation;
    shapeDyQ = (shapeTriDuQ * Jdudx(1,2) + shapeTriDvQ * Jdudx(2,2)) * orientation;

    % Source terms
    rhsQ = mySourceVolume(xQ, yQ);

    % Elemental matrices/vectors
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
        shapePhyQ = shapeRefLinQ * orientation;
        shapeAuxQ = shapePhyQ;

        % Mass matrices (physical space)
        weightsQ = weightsLinQ .* Jdxdu;
        matM_IIel = shapePhyQ' * (weightsQ .* shapePhyQ);
        matM_IGel = shapePhyQ' * (weightsQ .* shapeAuxQ);
        matM_GIel = shapeAuxQ' * (weightsQ .* shapePhyQ);
        matM_GGel = shapeAuxQ' * (weightsQ .* shapeAuxQ);

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
        switch FLUX

            case 'UPW'

                matIIel(dofLocP,dofLocP) = matIIel(dofLocP,dofLocP) + 1/eta(tri) * matM_IIel;
                matIIel(dofLocP,dofLocU) = matIIel(dofLocP,dofLocU) + nx * matM_IIel;
                matIIel(dofLocP,dofLocV) = matIIel(dofLocP,dofLocV) + ny * matM_IIel;
                matIGel = zeros(3*dofm.numDofPerTRI,dofm.numDofPerLIN);
                matIGel(dofLocP,:) = -1/eta(tri) * matM_IGel;
                matIGel(dofLocU,:) = nx * matM_IGel;
                matIGel(dofLocV,:) = ny * matM_IGel;

            case 'SYM'

                triNeigh = mesh.mapTriToTri(tri,fac);
                if (triNeigh>0)
                    etaF = sqrt(eta(tri)*eta(triNeigh));
                else
                    etaF = eta(tri);
                end

                matIIel(dofLocP,dofLocP) = matIIel(dofLocP,dofLocP) + 1/etaF * matM_IIel;
                matIIel(dofLocP,dofLocU) = matIIel(dofLocP,dofLocU) + nx * matM_IIel;
                matIIel(dofLocP,dofLocV) = matIIel(dofLocP,dofLocV) + ny * matM_IIel;
                matIGel = zeros(3*dofm.numDofPerTRI,dofm.numDofPerLIN);
                matIGel(dofLocP,:) = -1/etaF * matM_IGel;
                matIGel(dofLocU,:) = nx * matM_IGel;
                matIGel(dofLocV,:) = ny * matM_IGel;

            otherwise
                error('BAD FLUX.');
        end

        % -----------------------------------------------------------------
        % Auxiliary equations
        % -----------------------------------------------------------------

        % Elemental matrices/vectors
        matGGel = zeros(dofm.numDofPerLIN, dofm.numDofPerLIN);
        matGIel = zeros(dofm.numDofPerLIN, 3*dofm.numDofPerLIN);
        rhsGel  = zeros(dofm.numDofPerLIN, 1);

        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);

        if (triNeigh > 0)

            % Elemental matrices (interface condition)
            switch FLUX
                case 'UPW'
                    matGGel = 0.5 * matM_GGel;
                    matGIel = - eta(tri)*eta(triNeigh) / (eta(tri)+eta(triNeigh)) * [1/eta(tri) * matM_GIel, nx*matM_GIel, ny*matM_GIel];
                case 'SYM'
                    matGGel = 0.5 * matM_GGel;
                    etaF = sqrt(eta(tri)*eta(triNeigh));
                    matGIel = - etaF/2 * [1/etaF * matM_GIel, nx*matM_GIel, ny*matM_GIel];
                otherwise
                    error('BAD FLUX.');
            end

        else

            % Source terms
            [solP, ~, ~, solU, solV] = mySourceSurface(xQ, yQ);
            rhsPel = shapeAuxQ' * (weightsQ .* solP);
            rhsNUel = shapeAuxQ' * (weightsQ .* (nx.*solU + ny.*solV));

            % Type of BC
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));

            % Elemental matrices and RHS vectors (boundary conditions)
            switch FLUX
                case {'UPW','SYM'}
                    matGGel = matM_GGel;
                    switch BC
                        case 'DIR0'
                        case 'DIR'
                            rhsGel = rhsPel;
                        case 'NEU0'
                            matGIel = - [matM_GIel, eta(tri)*nx*matM_GIel, eta(tri)*ny*matM_GIel];
                        case 'NEU'
                            matGIel = - [matM_GIel, eta(tri)*nx*matM_GIel, eta(tri)*ny*matM_GIel];
                            rhsGel = -eta(tri)*rhsNUel;
                        case 'ABC'
                            matGIel = -1/2 * [matM_GIel, eta(tri)*nx*matM_GIel, eta(tri)*ny*matM_GIel];
                        case 'ROB'
                            matGIel = -1/2 * [matM_GIel, eta(tri)*nx*matM_GIel, eta(tri)*ny*matM_GIel];
                            rhsGel = (rhsPel - eta(tri)*rhsNUel) / 2;
                        otherwise
                            error('BAD BOUNDARY CONDITION.');
                    end
                otherwise
                    error('BAD FLUX.');
            end
        end

        % -----------------------------------------------------------------
        % Matrix assembling
        % -----------------------------------------------------------------

        % Global ID for edge unknowns
        edgGlo = abs(mesh.mapTriToEdg(tri,fac));
        dofGloG = dofm.locToGloLIN(edgGlo,:);
        dofLocG = 1:dofm.numDofPerLIN;
        if(mesh.mapTriToEdg(tri,fac) < 0)
            dofLocG(1) = 2;
            dofLocG(2) = 1;
        end

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

matGGvInv = zeros(numDofLIN, dofm.numDofPerLIN);
for lin=1:mesh.numEdg
    idLIN = (lin-1)*dofm.numDofPerLIN + (1:dofm.numDofPerLIN);
    matGGvInv(idLIN,:) = inv(matGGv(idLIN,:));
end

matII    = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG    = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofLIN);
matGI    = sparse(matGIx, matGIy, matGIv, numDofLIN, 3*numDofTRI);
matGG    = sparse(matGGx, matGGy, matGGv, numDofLIN, numDofLIN);
matIIinv = sparse(matIIx, matIIy, matIIvInv, 3*numDofTRI, 3*numDofTRI);
matGGinv = sparse(matGGx, matGGy, matGGvInv, numDofLIN, numDofLIN);

% -------------------------------------------------------------------------
% Source point
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

% Preconditioning
sysA.matP = matGG;
sysA.matPinv = matGGinv;

% Compute solution
solG = sysA.matS\sysA.rhsS;
solI = matIIinv*(rhsI-matIG*solG);

end
