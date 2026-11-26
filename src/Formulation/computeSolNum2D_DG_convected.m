function [solA, sysA] = computeSolNum2D_DG_convected(mesh, dofm)

global edgTagToBC
global omega rho c v0
global pntSouTag pntSouVal

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

TOT = 0;
SouEl = 0;

for tri=1:mesh.numTri

    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uTriQ, vTriQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5; % [ dx/du dx/dv ; dy/du dy/dv ]
    Jdudx = inv(Jdxdu); % [ du/dx du/dy ; dv/dx dv/dy ]
    detJdxdu = abs(det(Jdxdu));

    if(~isempty(pntSouTag))
        vertSou = mesh.mapPntToVer(mesh.tagPntFile == pntSouTag);
        for pos = 1:3
            if(mesh.mapTriToVer(tri,pos) == vertSou)
                TOT = TOT+1;
                SouEl(1,TOT) = tri;
                SouEl(2,TOT) = pos;
            end
        end
    end

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

    matMv1(dof,:) = -1i*omega/(rho*c^2)*matMel-1/(rho*c^2)*matDXel*v0(1)-1/(rho*c^2)*matDYel*v0(2);
    matMv2(dof,:) = -1i*omega*rho*matMel-rho*(matDXel*v0(1)+matDYel*v0(2));
    rhsP(dof) = -1/(1i*omega*rho) * vecRHSel;

end

matM1 = sparse(matXv,matYv,matMv1);
matM2 = sparse(matXv,matYv,matMv2);
matDX = sparse(matXv,matYv,matDXv);% Differentiation matrix (x)
matDY = sparse(matXv,matYv,matDYv);% Differentiation matrix (y)

matA = [
    matM1            -matDX                        -matDY                      ;
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
        Jdxdu = norm(V2-V1) * 0.5;% [ dx/du ]

        % Solution function
        %[solQ, solDxQ, solDyQ] = mySourceSurface(xQ, yQ);
        [solQ, ~, ~, solVxQ, solVyQ] = mySourceSurface(xQ,yQ);

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
        rhsUel = shapeOrQ' * weights * solVxQ * Jdxdu;
        rhsVel = shapeOrQ' * weights * solVyQ * Jdxdu;

        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);

        % Exterior tangent
        tx = -ny;
        ty = nx;

        v0n = v0(1)*nx+v0(2)*ny;

        if (v0n >= 0)
            gamma = 1;
        else
            gamma = 0;
        end

        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);

        if (triNeigh > 0)

            % Get global ID for exterior unknowns
            dofExt = dofLocTriNeigh(facNeigh,:);
            idExtP = dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = idExtP + numDofTRI;
            idExtV = idExtU + numDofTRI;

            matA(idIntP,idIntP) = matA(idIntP,idIntP) + (v0n+c)/(2*rho*c^2)                     * matMel;
            matA(idIntP,idIntU) = matA(idIntP,idIntU) + (v0n+c)/(2*c)                   *nx     * matMel;
            matA(idIntP,idIntV) = matA(idIntP,idIntV) + (v0n+c)/(2*c)                   *ny     * matMel;
            matA(idIntP,idExtP) = matA(idIntP,idExtP) + (v0n-c)/(2*rho*c^2)                     * matMel;
            matA(idIntP,idExtU) = matA(idIntP,idExtU) + (-v0n+c)/(2*c)                  *nx     * matMel;
            matA(idIntP,idExtV) = matA(idIntP,idExtV) + (-v0n+c)/(2*c)                  *ny     * matMel;

            matA(idIntU,idIntP) = matA(idIntU,idIntP) + (v0n+c)/(2*c)                   *nx     * matMel;
            matA(idIntU,idIntU) = matA(idIntU,idIntU) + (v0n+c)*rho/2                   *nx*nx  * matMel;
            matA(idIntU,idIntU) = matA(idIntU,idIntU) + v0n*rho            *gamma       *tx*tx  * matMel;
            matA(idIntU,idIntV) = matA(idIntU,idIntV) + (v0n+c)*rho/2                   *nx*ny  * matMel;
            matA(idIntU,idIntV) = matA(idIntU,idIntV) + v0n*rho            *gamma       *tx*ty  * matMel;
            matA(idIntU,idExtP) = matA(idIntU,idExtP) + (-v0n+c)/(2*c)                  *nx     * matMel;
            matA(idIntU,idExtU) = matA(idIntU,idExtU) + (v0n-c)*rho/2                   *nx*nx  * matMel;
            matA(idIntU,idExtU) = matA(idIntU,idExtU) + v0n*rho            *(1-gamma)   *tx*tx  * matMel;
            matA(idIntU,idExtV) = matA(idIntU,idExtV) + (v0n-c)*rho/2                   *nx*ny  * matMel;
            matA(idIntU,idExtV) = matA(idIntU,idExtV) + v0n*rho            *(1-gamma)   *tx*ty  * matMel;

            matA(idIntV,idIntP) = matA(idIntV,idIntP) + (v0n+c)/(2*c)                   *ny     * matMel;
            matA(idIntV,idIntU) = matA(idIntV,idIntU) + (v0n+c)*rho/2                   *nx*ny  * matMel;
            matA(idIntV,idIntU) = matA(idIntV,idIntU) + v0n*rho            *gamma       *tx*ty  * matMel;
            matA(idIntV,idIntV) = matA(idIntV,idIntV) + (v0n+c)*rho/2                   *ny*ny  * matMel;
            matA(idIntV,idIntV) = matA(idIntV,idIntV) + v0n*rho            *gamma       *ty*ty  * matMel;
            matA(idIntV,idExtP) = matA(idIntV,idExtP) + (-v0n+c)/(2*c)                  *ny     * matMel;
            matA(idIntV,idExtU) = matA(idIntV,idExtU) + (v0n-c)*rho/2                   *nx*ny  * matMel;
            matA(idIntV,idExtU) = matA(idIntV,idExtU) + v0n*rho            *(1-gamma)   *tx*ty  * matMel;
            matA(idIntV,idExtV) = matA(idIntV,idExtV) + (v0n-c)*rho/2                   *ny*ny  * matMel;
            matA(idIntV,idExtV) = matA(idIntV,idExtV) + v0n*rho            *(1-gamma)   *ty*ty  * matMel;

        else

            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));

            switch BC
                case 'DIR'

                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + 1/(rho*c) * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + nx        * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + ny        * matMel;

                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + v0n/c         *nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + v0n*rho    *nx*nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + v0n*rho    *tx*tx * matMel * gamma;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + v0n*rho    *ny*nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + v0n*rho    *ty*tx * matMel * gamma;

                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + v0n/c         *ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + v0n*rho    *nx*ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + v0n*rho    *tx*ty * matMel * gamma;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + v0n*rho    *ny*ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + v0n*rho    *ty*ty * matMel * gamma;

                    gp = rhsPel;
                    gux = rhsUel;
                    guy = rhsVel;
                    sR = gp;
                    rhsA(idIntP) = rhsA(idIntP) + sR * (c-v0n) / (rho*c^2);
                    rhsA(idIntU) = rhsA(idIntU) + sR * (-c+v0n) / c   * nx;
                    rhsA(idIntU) = rhsA(idIntU) - rho*v0n*(gux*tx+guy*ty)*tx * (1-gamma);
                    rhsA(idIntV) = rhsA(idIntV) + sR * (-c+v0n) / c   * ny;
                    rhsA(idIntV) = rhsA(idIntV) - rho*v0n*(gux*tx+guy*ty)*ty * (1-gamma);

                case 'NEU'

                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + v0n/(rho*c^2)     * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + v0n/c         *nx * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + v0n/c         *ny * matMel;

                    matA(idIntU,idIntP) = matA(idIntU,idIntP) +                nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + rho*c      *nx*nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + v0n*rho    *tx*tx * matMel * gamma;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + rho*c      *ny*nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + v0n*rho    *ty*tx * matMel * gamma;

                    matA(idIntV,idIntP) = matA(idIntV,idIntP) +                ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + rho*c      *nx*ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + v0n*rho    *tx*ty * matMel * gamma;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + rho*c      *ny*ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + v0n*rho    *ty*ty * matMel * gamma;

                    gux = rhsUel;
                    guy = rhsVel;
                    sN = gux*nx+guy*ny;
                    rhsA(idIntP) = rhsA(idIntP) + sN * (v0n-c) / c;
                    rhsA(idIntU) = rhsA(idIntU) + rho * sN * (c-v0n) * nx;
                    rhsA(idIntU) = rhsA(idIntU) - rho*v0n*(gux*tx+guy*ty)*tx * (1-gamma);
                    rhsA(idIntV) = rhsA(idIntV) + rho * sN * (c-v0n) * ny;
                    rhsA(idIntV) = rhsA(idIntV) - rho*v0n*(gux*tx+guy*ty)*ty * (1-gamma);

                case 'NEU0'

                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + v0n/(rho*c^2)     * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + v0n/c         *nx * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + v0n/c         *ny * matMel;

                    matA(idIntU,idIntP) = matA(idIntU,idIntP) +                nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + rho*c      *nx*nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + v0n*rho    *tx*tx * matMel * gamma;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + rho*c      *ny*nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + v0n*rho    *ty*tx * matMel * gamma;

                    matA(idIntV,idIntP) = matA(idIntV,idIntP) +                ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + rho*c      *nx*ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + v0n*rho    *tx*ty * matMel * gamma;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + rho*c      *ny*ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + v0n*rho    *ty*ty * matMel * gamma;

                    gux = rhsUel;
                    guy = rhsVel;
                    sN = 0*(gux*nx+guy*ny);
                    rhsA(idIntP) = rhsA(idIntP) + sN * (v0n-c) / c;
                    rhsA(idIntU) = rhsA(idIntU) + rho * sN * (c-v0n) * nx;
                    rhsA(idIntU) = rhsA(idIntU) - rho*v0n*(gux*tx+guy*ty)*tx * (1-gamma);
                    rhsA(idIntV) = rhsA(idIntV) + rho * sN * (c-v0n) * ny;
                    rhsA(idIntV) = rhsA(idIntV) - rho*v0n*(gux*tx+guy*ty)*ty * (1-gamma);

                case 'ABC'

                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + (v0n+c)/(2*rho*c^2)     * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + (v0n+c)/(2*c)       *nx * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + (v0n+c)/(2*c)       *ny * matMel;

                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + (v0n+c)/(2*c)       *nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + rho*(v0n+c)/2    *nx*nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + v0n*rho          *tx*tx * matMel * gamma;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + rho*(v0n+c)/2    *ny*nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + v0n*rho          *ty*tx * matMel * gamma;

                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + (v0n+c)/(2*c)       *ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + rho*(v0n+c)/2    *nx*ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + v0n*rho          *tx*ty * matMel * gamma;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + rho*(v0n+c)/2    *ny*ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + v0n*rho          *ty*ty * matMel * gamma;

                    gp = rhsPel*0;
                    gux = rhsUel*0;
                    guy = rhsVel*0;
                    sR = 0*(gp-rho*c*(gux*nx+guy*ny));
                    rhsA(idIntP) = rhsA(idIntP) + sR * (c-v0n) / (2*rho*c^2);
                    rhsA(idIntU) = rhsA(idIntU) + sR * (v0n-c)/(2*c) * nx;
                    rhsA(idIntU) = rhsA(idIntU) - 0*rho*v0n*(gux*tx+guy*ty)*tx * (1-gamma);
                    rhsA(idIntV) = rhsA(idIntV) + sR * (v0n-c)/(2*c) * ny;
                    rhsA(idIntV) = rhsA(idIntV) - 0*rho*v0n*(gux*tx+guy*ty)*ty * (1-gamma);

                case 'ROB'

                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + (v0n+c)/(2*rho*c^2)     * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + (v0n+c)/(2*c)       *nx * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + (v0n+c)/(2*c)       *ny * matMel;

                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + (v0n+c)/(2*c)       *nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + rho*(v0n+c)/2    *nx*nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + v0n*rho          *tx*tx * matMel * gamma;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + rho*(v0n+c)/2    *ny*nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + v0n*rho          *ty*tx * matMel * gamma;

                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + (v0n+c)/(2*c)       *ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + rho*(v0n+c)/2    *nx*ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + v0n*rho          *tx*ty * matMel * gamma;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + rho*(v0n+c)/2    *ny*ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + v0n*rho          *ty*ty * matMel * gamma;

                    gp = rhsPel;
                    gux = rhsUel;
                    guy = rhsVel;
                    sR = gp-rho*c*(gux*nx+guy*ny);
                    rhsA(idIntP) = rhsA(idIntP) + sR * (c-v0n) / (2*rho*c^2);
                    rhsA(idIntU) = rhsA(idIntU) + sR * (v0n-c)/(2*c) * nx;
                    rhsA(idIntU) = rhsA(idIntU) - rho*v0n*(gux*tx+guy*ty)*tx * (1-gamma);
                    rhsA(idIntV) = rhsA(idIntV) + sR * (v0n-c)/(2*c) * ny;
                    rhsA(idIntV) = rhsA(idIntV) - rho*v0n*(gux*tx+guy*ty)*ty * (1-gamma);

                case 'GC'

                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + (v0n+c)/(2*rho*c^2)     * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + (v0n+c)/(2*c)       *nx * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + (v0n+c)/(2*c)       *ny * matMel;

                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + (v0n+c)/(2*c)       *nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + rho*(v0n+c)/2    *nx*nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + v0n*rho          *tx*tx * matMel * gamma;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + rho*(v0n+c)/2    *ny*nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + v0n*rho          *ty*tx * matMel * gamma;

                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + (v0n+c)/(2*c)       *ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + rho*(v0n+c)/2    *nx*ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + v0n*rho          *tx*ty * matMel * gamma;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + rho*(v0n+c)/2    *ny*ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + v0n*rho          *ty*ty * matMel * gamma;

                    gp = rhsPel;
                    gux = rhsUel;
                    guy = rhsVel;
                    rhsA(idIntP) = rhsA(idIntP) + (c-v0n) / (2*rho*c^2) * gp + (v0n-c) / (2*c) * (gux*nx+guy*ny);
                    rhsA(idIntU) = rhsA(idIntU) + ((v0n-c)/(2*c)*gp + rho*(c-v0n)/2*(gux*nx+guy*ny))*nx;
                    rhsA(idIntU) = rhsA(idIntU) - rho*v0n*(gux*tx+guy*ty)*tx * (1-gamma);
                    rhsA(idIntV) = rhsA(idIntV) + ((v0n-c)/(2*c)*gp + rho*(c-v0n)/2*(gux*nx+guy*ny))*ny;
                    rhsA(idIntV) = rhsA(idIntV) - rho*v0n*(gux*tx+guy*ty)*ty * (1-gamma);

                otherwise
                    error('BAD BOUNDARY CONDITION.');
            end
        end
    end
end

if(~isempty(pntSouTag))
    for ind=1:TOT
        dofSou = dofm.numDofPerTRI * (SouEl(1,ind)-1) + SouEl(2,ind);
        rhsA(dofSou) = rhsA(dofSou) +  pntSouVal / TOT;
    end
end

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

% Full system
sysA.matA = matA;
sysA.rhsA = rhsA;

% Preconditionning - NO PRECONDITIONNING TECHNIQUE CODED YET FOR DG.
sysA.matP = 1;
sysA.matPinv = 1;

% Reduced system
sysA.matS = matA;
sysA.rhsS = rhsA;

% Physical system
sysA.matPhy = matA;
sysA.rhsPhy = rhsA;

% Compute solution
solA = matA\rhsA;
% solA = rhsA * 0;

end