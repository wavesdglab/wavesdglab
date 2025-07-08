% Contributors: Simone Pescuma, Axel Modave

function [solA, sysA] = computeSolNum2D_DG_heterogeneous(mesh, dofm, PREC)

global k eta edgTagToBC


numDofTRI = dofm.numDofTRI;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);
weightsTriQ = sparse(1:size(weightsTriQ,1), 1:size(weightsTriQ,1), weightsTriQ);

% Shape functions (f, dfdu, dfdv)
shapeQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeDuQ, shapeDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

% Global matrices
matXv  = zeros(numDofTRI, numDofPerTRI);
matYv  = zeros(numDofTRI, numDofPerTRI);
matMv1  = zeros(numDofTRI, numDofPerTRI);
matMv2  = zeros(numDofTRI, numDofPerTRI);
matDXv = zeros(numDofTRI, numDofPerTRI);
matDYv = zeros(numDofTRI, numDofPerTRI);
rhsP   = zeros(numDofTRI, 1);

for tri=1:mesh.numTri

    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uTriQ, vTriQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    Jdudx = inv(Jdxdu);                 % [ du/dx du/dy ; dv/dx dv/dy ]
    detJdxdu = abs(det(Jdxdu));

    % Orientation
    orientation = ones(numDofPerTRI,1);
    if(ver(1) > ver(2))
        orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(ver(2) > ver(3))
        orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(ver(3) > ver(1))
        orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    orientation = sparse(1:numDofPerTRI, 1:numDofPerTRI, orientation);

    % Shape functions (f, dfdx, dfdy) with orientation
    shapeOrQ = shapeQ * orientation;
    shapeDxQ = (shapeDuQ * Jdudx(1,1) + shapeDvQ * Jdudx(2,1)) * orientation;
    shapeDyQ = (shapeDuQ * Jdudx(1,2) + shapeDvQ * Jdudx(2,2)) * orientation;

    % RHS function
    rhsQ = mySourceVolume(xQ, yQ);

    % Elemental matrices
    matMel = shapeOrQ' * weightsTriQ * shapeOrQ * detJdxdu;
    matDXel = shapeDxQ' * weightsTriQ * shapeOrQ * detJdxdu;
    matDYel = shapeDyQ' * weightsTriQ * shapeOrQ * detJdxdu;
    vecRHSel = shapeOrQ' * weightsTriQ * rhsQ * detJdxdu;

    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matXv(dof,:) = dof'*ones(1,size(dof,2));
    matYv(dof,:) = ones(size(dof,2),1)*dof;
    matDXv(dof,:) = matDXel;
    matDYv(dof,:) = matDYel;

    matMv1(dof,:) = -1i*k(tri)/eta(tri)  * matMel;
    matMv2(dof,:) = -1i*k(tri)*eta(tri)  * matMel;
    rhsP(dof) = -1/(1i*k(tri)*eta(tri)) * vecRHSel;

end

matM1 = sparse(matXv,matYv,matMv1);
matM2 = sparse(matXv,matYv,matMv2);
matDX = sparse(matXv,matYv,matDXv);  % Differentiation matrix (x)
matDY = sparse(matXv,matYv,matDYv);  % Differentiation matrix (y)

matA = [
    matM1           -matDX                        -matDY                      ;
    -matDX           matM2                         sparse(numDofTRI,numDofTRI) ;
    -matDY           sparse(numDofTRI,numDofTRI)   matM2                        ];

rhsA = [
    rhsP ;
    zeros(numDofTRI,1) ;
    zeros(numDofTRI,1) ];

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

% Quadrature
degreeQ = 4*dofm.degree;
[uQ, weights] = quadratureGaussLIN(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);

% Shape functions
shapeQ = functionsShapeLIN(uQ, dofm.degree);

dofLocTri = zeros(3,dofm.numDofPerLIN);
dofLocTri(1,:) = [1 2 dofm.locEdg(1,:)];
dofLocTri(2,:) = [2 3 dofm.locEdg(2,:)];
dofLocTri(3,:) = [3 1 dofm.locEdg(3,:)];

dofLocTriNeigh = zeros(3,dofm.numDofPerLIN);
dofLocTriNeigh(1,:) = [2 1 dofm.locEdg(1,:)];
dofLocTriNeigh(2,:) = [3 2 dofm.locEdg(2,:)];
dofLocTriNeigh(3,:) = [1 3 dofm.locEdg(3,:)];

for tri=1:mesh.numTri

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

        % Global ID for interior unknowns
        dofInt = dofLocTri(fac,:);
        idIntP = 0*numDofTRI + dofm.locToGloTRI(tri,dofInt);
        idIntU = 1*numDofTRI + dofm.locToGloTRI(tri,dofInt);
        idIntV = 2*numDofTRI + dofm.locToGloTRI(tri,dofInt);

        % Mapping
        V1 = mesh.coord(n1(fac),:);
        V2 = mesh.coord(n2(fac),:);
        [xQ, yQ] = locToGloLIN(uQ, V1, V2);
        Jdxdu = norm(V2-V1) * 0.5;  % [ dx/du ]

        % Solution function
        [solQ, solDxQ, solDyQ] = mySourceSurface(xQ, yQ);

        % Orientation
        orientation = ones(dofm.numDofPerLIN,1);
        if(n1(fac) > n2(fac))
            orientation(3:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        orientation = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation);

        % Shape function
        shapeOrQ = shapeQ * orientation;

        % Elemental matrix
        matMel = shapeOrQ' * weights * shapeOrQ * Jdxdu;
        rhsPel = shapeOrQ' * weights * solQ * Jdxdu;
        rhsUel = shapeOrQ' * weights * solDxQ * Jdxdu / (1i*k(tri)*eta(tri));
        rhsVel = shapeOrQ' * weights * solDyQ * Jdxdu / (1i*k(tri)*eta(tri));

        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);

        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);

        if (triNeigh > 0)

            % Get global ID for exterior unknowns
            dofExt = dofLocTriNeigh(facNeigh,:);
            idExtP = dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = idExtP + numDofTRI;
            idExtV = idExtU + numDofTRI;

            matA(idIntP,idIntP) = matA(idIntP,idIntP) + 1/(eta(tri)+eta(triNeigh))                     * matMel;
            matA(idIntP,idIntU) = matA(idIntP,idIntU) + eta(tri)/(eta(tri)+eta(triNeigh))           *nx     * matMel;
            matA(idIntP,idIntV) = matA(idIntP,idIntV) + eta(tri)/(eta(tri)+eta(triNeigh))           *ny     * matMel;
            matA(idIntP,idExtP) = matA(idIntP,idExtP) - 1/(eta(tri)+eta(triNeigh))                     * matMel;
            matA(idIntP,idExtU) = matA(idIntP,idExtU) + eta(triNeigh)/(eta(tri)+eta(triNeigh))      *nx     * matMel;
            matA(idIntP,idExtV) = matA(idIntP,idExtV) + eta(triNeigh)/(eta(tri)+eta(triNeigh))      *ny     * matMel;

            matA(idIntU,idIntP) = matA(idIntU,idIntP) + eta(triNeigh)/(eta(tri)+eta(triNeigh))      *nx     * matMel;
            matA(idIntU,idIntU) = matA(idIntU,idIntU) + eta(tri)*eta(triNeigh)/(eta(tri)+eta(triNeigh))  *nx*nx  * matMel;
            matA(idIntU,idIntV) = matA(idIntU,idIntV) + eta(tri)*eta(triNeigh)/(eta(tri)+eta(triNeigh))  *nx*ny  * matMel;
            matA(idIntU,idExtP) = matA(idIntU,idExtP) + eta(tri)/(eta(tri)+eta(triNeigh))           *nx     * matMel;
            matA(idIntU,idExtU) = matA(idIntU,idExtU) - eta(tri)*eta(triNeigh)/(eta(tri)+eta(triNeigh))  *nx*nx  * matMel;
            matA(idIntU,idExtV) = matA(idIntU,idExtV) - eta(tri)*eta(triNeigh)/(eta(tri)+eta(triNeigh))  *nx*ny  * matMel;

            matA(idIntV,idIntP) = matA(idIntV,idIntP) + eta(triNeigh)/(eta(tri)+eta(triNeigh))      *ny     * matMel;
            matA(idIntV,idIntU) = matA(idIntV,idIntU) + eta(tri)*eta(triNeigh)/(eta(tri)+eta(triNeigh))  *nx*ny  * matMel;
            matA(idIntV,idIntV) = matA(idIntV,idIntV) + eta(tri)*eta(triNeigh)/(eta(tri)+eta(triNeigh))  *ny*ny  * matMel;
            matA(idIntV,idExtP) = matA(idIntV,idExtP) + eta(tri)/(eta(tri)+eta(triNeigh))           *ny     * matMel;
            matA(idIntV,idExtU) = matA(idIntV,idExtU) - eta(tri)*eta(triNeigh)/(eta(tri)+eta(triNeigh))  *nx*ny  * matMel;
            matA(idIntV,idExtV) = matA(idIntV,idExtV) - eta(tri)*eta(triNeigh)/(eta(tri)+eta(triNeigh))  *ny*ny  * matMel;

        else

            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));

            switch BC
                case 'DIR'

                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + 1/eta(tri) * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + nx        * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + ny        * matMel;

                    gp = rhsPel;
                    rhsA(idIntP) = rhsA(idIntP) + gp / eta(tri);
                    rhsA(idIntU) = rhsA(idIntU) - gp * nx;
                    rhsA(idIntV) = rhsA(idIntV) - gp * ny;

                case 'NEU'

                    matA(idIntU,idIntP) = matA(idIntU,idIntP) +                  nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + eta(tri) * nx       * nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + eta(tri) * ny       * nx * matMel;

                    matA(idIntV,idIntP) = matA(idIntV,idIntP) +                  ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + eta(tri) * nx       * ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + eta(tri) * ny       * ny * matMel;

                    gnu = nx*rhsUel + ny*rhsVel;
                    rhsA(idIntP) = rhsA(idIntP) - gnu;
                    rhsA(idIntU) = rhsA(idIntU) + eta(tri) * gnu * nx;
                    rhsA(idIntV) = rhsA(idIntV) + eta(tri) * gnu * ny;

                case 'ABC'

                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5 / eta(tri)           * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5       * nx      * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + 0.5       * ny      * matMel;

                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5            * nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5 * eta(tri) * nx * nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5 * eta(tri) * ny * nx * matMel;

                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + 0.5            * ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5 * eta(tri) * nx * ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5 * eta(tri) * ny * ny * matMel;

                    gchar = rhsPel - etaNeigh * (nx*rhsUel + ny*rhsVel);
                    rhsA(idIntP) = rhsA(idIntP) + gchar * 0.5 / eta(tri);
                    rhsA(idIntU) = rhsA(idIntU) - gchar * 0.5 * nx;
                    rhsA(idIntV) = rhsA(idIntV) - gchar * 0.5 * ny;

                case 'ROB'

                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5 / eta(tri)           * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5       * nx      * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + 0.5       * ny      * matMel;

                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5            * nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5 * eta(tri) * nx * nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5 * eta(tri) * ny * nx * matMel;

                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + 0.5            * ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5 * eta(tri) * nx * ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5 * eta(tri) * ny * ny * matMel;

                    gchar = rhsPel - etaNeigh * (nx*rhsUel + ny*rhsVel);
                    rhsA(idIntP) = rhsA(idIntP) + gchar * 0.5 / eta(tri);
                    rhsA(idIntU) = rhsA(idIntU) - gchar * 0.5 * nx;
                    rhsA(idIntV) = rhsA(idIntV) - gchar * 0.5 * ny;

                otherwise
                    error('BAD BOUNDARY CONDITION.');
            end
        end
    end
end

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

% Full system
sysA.matA = matA;
sysA.rhsA = rhsA;

% Preconditionning
if (PREC == 1)
    warning('NO PRECONDITIONNING TECHNIQUE CODED YET FOR DG.')
    sysA.matP = matP;
    sysA.matPinv = matPinv;
else
    sysA.matP = 1;
    sysA.matPinv = 1;
end

% Compute solution
solA = matA\rhsA;

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
