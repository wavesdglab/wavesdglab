% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solA, sysA] = computeSolNum2D_DG(mesh, dofm, tau, theta, PREC)

global k edgTagToBC
global LdomX LdomY LpmlX LpmlY Rdom Rpml

numDofTRI = dofm.numDofTRI;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);

% Shape functions and derivatives (reference space)
shapeLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeTriDuQ, shapeTriDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

% Global matrices
matXv  = zeros(numDofTRI, numDofPerTRI);
matYv  = zeros(numDofTRI, numDofPerTRI);
matMv  = zeros(numDofTRI, numDofPerTRI);
matMPv = zeros(numDofTRI, numDofPerTRI);
matMXXv = zeros(numDofTRI, numDofPerTRI);
matMXYv = zeros(numDofTRI, numDofPerTRI);
matMYYv = zeros(numDofTRI, numDofPerTRI);
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
    shapeOrQ = shapeTriQ * orientation;
    shapeDxQ = (shapeTriDuQ * Jdudx(1,1) + shapeTriDvQ * Jdudx(2,1)) * orientation;
    shapeDyQ = (shapeTriDuQ * Jdudx(1,2) + shapeTriDvQ * Jdudx(2,2)) * orientation;

    % RHS function
    [~, ~, ~, rhsQ] = mySol(xQ, yQ);

    % Elemental matrices
    weightsQ = weightsTriQ .* detJdxdu;
    matMel = transpose(shapeOrQ) * (weightsQ .* shapeOrQ);
    matMPel = matMel;
    matMXXel = matMel;
    matMXYel = zeros(size(matMel));
    matMYYel = matMel;
    matDXel = transpose(shapeDxQ) * (weightsQ .* shapeOrQ);
    matDYel = transpose(shapeDyQ) * (weightsQ .* shapeOrQ);
    vecRHSel = transpose(shapeOrQ) * (weightsQ .* rhsQ);

    % PML stretching (rectangular PML)
    if(~isempty(LdomX) && ~isempty(LdomY))
        if ((mean(abs(xQ)) >= LdomX) || (mean(abs(yQ)) >= LdomY))
            sigmaPmlX = 1./(LdomX+LpmlX-abs(xQ));
            sigmaPmlY = 1./(LdomY+LpmlY-abs(yQ));
            gammaPmlX = ones(size(xQ)) - sigmaPmlX/(1i*k);
            gammaPmlY = ones(size(yQ)) - sigmaPmlY/(1i*k);
            coefPml = gammaPmlX.*gammaPmlY;
            tensPmlInvXX = gammaPmlX./gammaPmlY;
            tensPmlInvYY = gammaPmlY./gammaPmlX;
            matMPel = transpose(shapeOrQ) * (weightsQ .* coefPml .* shapeOrQ);
            matMXXel = transpose(shapeOrQ) * (weightsQ .* tensPmlInvXX .* shapeOrQ);
            matMYYel = transpose(shapeOrQ) * (weightsQ .* tensPmlInvYY .* shapeOrQ);
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
            matMPel = transpose(shapeOrQ) * (weightsQ .* coefPml .* shapeOrQ);
            matMXXel = transpose(shapeOrQ) * (weightsQ .* tensPmlInvXX .* shapeOrQ);
            matMXYel = transpose(shapeOrQ) * (weightsQ .* tensPmlInvXY .* shapeOrQ);
            matMYYel = transpose(shapeOrQ) * (weightsQ .* tensPmlInvYY .* shapeOrQ);
        end
    end

    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matXv(dof,:) = dof'*ones(1,size(dof,2));
    matYv(dof,:) = ones(size(dof,2),1)*dof;
    matMv(dof,:) = matMel;
    matMPv(dof,:) = matMPel;
    matMXXv(dof,:) = matMXXel;
    matMXYv(dof,:) = matMXYel;
    matMYYv(dof,:) = matMYYel;
    matDXv(dof,:) = matDXel;
    matDYv(dof,:) = matDYel;
    rhsP(dof) = vecRHSel;

end

matM  = sparse(matXv,matYv,matMv);     % Mass matrix
matMP = sparse(matXv,matYv,matMPv);    % Mass matrix
matMXX = sparse(matXv,matYv,matMXXv);  % Mass matrix
matMXY = sparse(matXv,matYv,matMXYv);  % Mass matrix
matMYY = sparse(matXv,matYv,matMYYv);  % Mass matrix
matDX = sparse(matXv,matYv,matDXv);    % Differentiation matrix (x)
matDY = sparse(matXv,matYv,matDYv);    % Differentiation matrix (y)

matP = [
    matM sparse(numDofTRI,numDofTRI) sparse(numDofTRI,numDofTRI) ;
    sparse(numDofTRI,numDofTRI) matM sparse(numDofTRI,numDofTRI) ;
    sparse(numDofTRI,numDofTRI) sparse(numDofTRI,numDofTRI) matM ];

matA = [
    -1i*k*matMP  -matDX        -matDY       ;
    -matDX       -1i*k*matMXX  -1i*k*matMXY ;
    -matDY       -1i*k*matMXY  -1i*k*matMYY ];

rhsA = [
    -1/(1i*k)*rhsP ;
    zeros(numDofTRI,1) ;
    zeros(numDofTRI,1) ];

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

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
        [xQ, yQ] = locToGloLIN(uLinQ, V1, V2);
        detJdxdu = norm(V2-V1) * 0.5;  % [ dx/du ]

        % Solution function
        [solQ, solDxQ, solDyQ, ~] = mySol(xQ, yQ);

        % Orientation
        orientation = ones(dofm.numDofPerLIN,1);
        if(n1(fac) > n2(fac))
            orientation(3:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        orientation = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation);

        % Shape function
        shapeOrQ = shapeLinQ * orientation;

        % Elemental matrix
        weightsQ = weightsLinQ .* detJdxdu;
        matMel = transpose(shapeOrQ) * (weightsQ .* shapeOrQ);
        rhsPel = transpose(shapeOrQ) * (weightsQ .* solQ);
        rhsUel = transpose(shapeOrQ) * (weightsQ .* solDxQ / (1i*k));
        rhsVel = transpose(shapeOrQ) * (weightsQ .* solDyQ / (1i*k));

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

            matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau*theta   * matMel;
            matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5*nx          * matMel;
            matA(idIntP,idIntV) = matA(idIntP,idIntV) + 0.5*ny          * matMel;
            matA(idIntP,idExtP) = matA(idIntP,idExtP) - 0.5*tau*theta   * matMel;
            matA(idIntP,idExtU) = matA(idIntP,idExtU) + 0.5*nx          * matMel;
            matA(idIntP,idExtV) = matA(idIntP,idExtV) + 0.5*ny          * matMel;

            matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5*nx                * matMel;
            matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau*theta * nx*nx * matMel;
            matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5/tau*theta * nx*ny * matMel;
            matA(idIntU,idExtP) = matA(idIntU,idExtP) + 0.5*nx                * matMel;
            matA(idIntU,idExtU) = matA(idIntU,idExtU) - 0.5/tau*theta * nx*nx * matMel;
            matA(idIntU,idExtV) = matA(idIntU,idExtV) - 0.5/tau*theta * nx*ny * matMel;

            matA(idIntV,idIntP) = matA(idIntV,idIntP) + 0.5*ny                * matMel;
            matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5/tau*theta * nx*ny * matMel;
            matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5/tau*theta * ny*ny * matMel;
            matA(idIntV,idExtP) = matA(idIntV,idExtP) + 0.5*ny                * matMel;
            matA(idIntV,idExtU) = matA(idIntV,idExtU) - 0.5/tau*theta * nx*ny * matMel;
            matA(idIntV,idExtV) = matA(idIntV,idExtV) - 0.5/tau*theta * ny*ny * matMel;

        else

            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = edgTagToBC(mesh.tagEdg(edgGlo));

            switch BC
                case {'DIR0','DIR'}

                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + tau*theta * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + nx        * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + ny        * matMel;

                    if strcmp(BC,'DIR')
                        gp = rhsPel;
                        rhsA(idIntP) = rhsA(idIntP) + gp * tau*theta;
                        rhsA(idIntU) = rhsA(idIntU) - gp * nx;
                        rhsA(idIntV) = rhsA(idIntV) - gp * ny;
                    end

                case {'NEU0','NEU'}

                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + 1            * nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + nx/tau*theta * nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + ny/tau*theta * nx * matMel;

                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + 1            * ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + nx/tau*theta * ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + ny/tau*theta * ny * matMel;

                    if strcmp(BC,'NEU')
                        gnu = nx*rhsUel + ny*rhsVel;
                        rhsA(idIntP) = rhsA(idIntP) - gnu;
                        rhsA(idIntU) = rhsA(idIntU) + gnu * nx/tau*theta;
                        rhsA(idIntV) = rhsA(idIntV) + gnu * ny/tau*theta;
                    end

                case {'ABC','ROB'}

                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + tau/(1+tau)         * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + 1/(1+tau) * nx      * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + 1/(1+tau) * ny      * matMel;

                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + tau/(1+tau)    * nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + 1/(1+tau) * nx * nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + 1/(1+tau) * ny * nx * matMel;

                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + tau/(1+tau)    * ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + 1/(1+tau) * nx * ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + 1/(1+tau) * ny * ny * matMel;

                    if strcmp(BC,'ROB')
                        gchar = rhsPel - (nx*rhsUel + ny*rhsVel);
                        rhsA(idIntP) = rhsA(idIntP) + gchar * tau/(1+tau);
                        rhsA(idIntU) = rhsA(idIntU) - gchar * 1/(1+tau) * nx;
                        rhsA(idIntV) = rhsA(idIntV) - gchar * 1/(1+tau) * ny;
                    end

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
    sysA.matP = matP;
else
    sysA.matP = 1;
end

% Compute solution
solA = matA\rhsA;

end