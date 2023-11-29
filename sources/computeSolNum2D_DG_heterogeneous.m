% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [solA, sysA] = computeSolNum2D_DG_heterogeneous(mesh, dofm, PREC)

global omega

numDofTRI = dofm.numDofTRI;
numDofPerTRI = dofm.numDofPerTRI;

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

% Quadrature
degreeQ = 2*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);

% Shape functions (f, dfdu, dfdv)
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);
[shapeDuQ, shapeDvQ] = functionsShapeDerTRI(uQ, vQ, dofm.degree);

% Global matrices
matXv  = zeros(numDofTRI, numDofPerTRI);
matYv  = zeros(numDofTRI, numDofPerTRI);
matMv  = zeros(numDofTRI, numDofPerTRI);
matMvInv = zeros(numDofTRI, numDofPerTRI);
matDXv = zeros(numDofTRI, numDofPerTRI);
matDYv = zeros(numDofTRI, numDofPerTRI);
rhsP   = zeros(numDofTRI, 1);

for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);  
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uQ, vQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    Jdudx = inv(Jdxdu);                 % [ du/dx du/dy ; dv/dx dv/dy ]
    detJdxdu = abs(det(Jdxdu));
    
    % Physical parameters on the element
    x_C = (V1(1,1)+V2(1,1)+V3(1,1))/3;
    y_C = (V1(1,2)+V2(1,2)+V3(1,2))/3;
    [~, ~, ~, ~, ~, ~, ~, c, eta] = mySol2D_heterogeneous(x_C,y_C);

    k = omega / c;

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
    [~, ~, ~, rhsQ, ~, ~, ~, ~, ~] = mySol2D_heterogeneous(xQ,yQ);
    
    % Elemental matrices
    matMel = shapeOrQ' * weights * shapeOrQ * detJdxdu;
    matDXel = shapeDxQ' * weights * shapeOrQ * detJdxdu;
    matDYel = shapeDyQ' * weights * shapeOrQ * detJdxdu;
    vecRHSel = shapeOrQ' * weights * rhsQ * detJdxdu;
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matXv(dof,:) = dof'*ones(1,size(dof,2));
    matYv(dof,:) = ones(size(dof,2),1)*dof;
    matMv(dof,:) = matMel;
    matMvInv(dof,:) = inv(matMel);
    matDXv(dof,:) = matDXel;
    matDYv(dof,:) = matDYel;
    rhsP(dof) = vecRHSel;
    
end

matM  = sparse(matXv,matYv,matMv);   % Mass matrix
matMinv  = sparse(matXv,matYv,matMvInv); % Mass matrix (inverse)
matDX = sparse(matXv,matYv,matDXv);  % Differentiation matrix (x)
matDY = sparse(matXv,matYv,matDYv);  % Differentiation matrix (y)

matP = [
    matM sparse(numDofTRI,numDofTRI) sparse(numDofTRI,numDofTRI) ;
    sparse(numDofTRI,numDofTRI) matM sparse(numDofTRI,numDofTRI) ;
    sparse(numDofTRI,numDofTRI) sparse(numDofTRI,numDofTRI) matM ];

matPinv = [
    matMinv sparse(numDofTRI,numDofTRI) sparse(numDofTRI,numDofTRI) ;
    sparse(numDofTRI,numDofTRI) matMinv sparse(numDofTRI,numDofTRI) ;
    sparse(numDofTRI,numDofTRI) sparse(numDofTRI,numDofTRI) matMinv ];

matA = [
    -1i*k/eta*matM   -matDX                        -matDY                      ;
    -matDX           -1i*k*eta*matM                sparse(numDofTRI,numDofTRI) ;
    -matDY           sparse(numDofTRI,numDofTRI)   -1i*k*eta*matM              ];

rhsA = [
    -1/(1i*k*eta)*rhsP ;
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
        [solQ, solDxQ, solDyQ, ~, ~, ~, ~, ~, ~] = mySol2D_heterogeneous(xQ,yQ);
        
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
        rhsUel = shapeOrQ' * weights * solDxQ * Jdxdu / (1i*k*eta);
        rhsVel = shapeOrQ' * weights * solDyQ * Jdxdu / (1i*k*eta);
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);

        if (triNeigh>0)
            verTri = mesh.mapTriToVer(triNeigh,:);
            V1 = mesh.coord(verTri(1),:);
            V2 = mesh.coord(verTri(2),:);
            V3 = mesh.coord(verTri(3),:);

            % Physical parameters on the neighboring element
            x_C = (V1(1,1)+V2(1,1)+V3(1,1))/3;
            y_C = (V1(1,2)+V2(1,2)+V3(1,2))/3;
            [~, ~, ~, ~, ~, ~, ~, ~, etaNeigh] = mySol2D_heterogeneous(x_C,y_C);
        else
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            BC = tagToBC(mesh.tagEdg(edgGlo));
            switch BC
                case {'DIR', 'NEU', 'ABC'}
                    etaNeigh = eta;
                otherwise
                    error('BAD BOUNDARY CONDITION.');
            end
        end

        if (triNeigh > 0)
            
            % Get global ID for exterior unknowns
            dofExt = dofLocTriNeigh(facNeigh,:);
            idExtP = dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = idExtP + numDofTRI;
            idExtV = idExtU + numDofTRI;
            
            matA(idIntP,idIntP) = matA(idIntP,idIntP) + 1/(eta+etaNeigh)                     * matMel;
            matA(idIntP,idIntU) = matA(idIntP,idIntU) + eta/(eta+etaNeigh)           *nx     * matMel;
            matA(idIntP,idIntV) = matA(idIntP,idIntV) + eta/(eta+etaNeigh)           *ny     * matMel;
            matA(idIntP,idExtP) = matA(idIntP,idExtP) - 1/(eta+etaNeigh)                     * matMel;
            matA(idIntP,idExtU) = matA(idIntP,idExtU) + etaNeigh/(eta+etaNeigh)      *nx     * matMel;
            matA(idIntP,idExtV) = matA(idIntP,idExtV) + etaNeigh/(eta+etaNeigh)      *ny     * matMel;
            
            matA(idIntU,idIntP) = matA(idIntU,idIntP) + etaNeigh/(eta+etaNeigh)      *nx     * matMel;
            matA(idIntU,idIntU) = matA(idIntU,idIntU) + eta*etaNeigh/(eta+etaNeigh)  *nx*nx  * matMel;
            matA(idIntU,idIntV) = matA(idIntU,idIntV) + eta*etaNeigh/(eta+etaNeigh)  *nx*ny  * matMel;
            matA(idIntU,idExtP) = matA(idIntU,idExtP) + eta/(eta+etaNeigh)           *nx     * matMel;
            matA(idIntU,idExtU) = matA(idIntU,idExtU) - eta*etaNeigh/(eta+etaNeigh)  *nx*nx  * matMel;
            matA(idIntU,idExtV) = matA(idIntU,idExtV) - eta*etaNeigh/(eta+etaNeigh)  *nx*ny  * matMel;
            
            matA(idIntV,idIntP) = matA(idIntV,idIntP) + etaNeigh/(eta+etaNeigh)      *ny     * matMel;
            matA(idIntV,idIntU) = matA(idIntV,idIntU) + eta*etaNeigh/(eta+etaNeigh)  *nx*ny  * matMel;
            matA(idIntV,idIntV) = matA(idIntV,idIntV) + eta*etaNeigh/(eta+etaNeigh)  *ny*ny  * matMel;
            matA(idIntV,idExtP) = matA(idIntV,idExtP) + eta/(eta+etaNeigh)           *ny     * matMel;
            matA(idIntV,idExtU) = matA(idIntV,idExtU) - eta*etaNeigh/(eta+etaNeigh)  *nx*ny  * matMel;
            matA(idIntV,idExtV) = matA(idIntV,idExtV) - eta*etaNeigh/(eta+etaNeigh)  *ny*ny  * matMel;
            
        else
            
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            switch tagToBC(mesh.tagEdg(edgGlo))
                case 'DIR'
                    
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + 1/eta     * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + nx        * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + ny        * matMel;
                    
                    gp = rhsPel;
                    rhsA(idIntP) = rhsA(idIntP) + gp / eta;
                    rhsA(idIntU) = rhsA(idIntU) - gp * nx;
                    rhsA(idIntV) = rhsA(idIntV) - gp * ny;
                    
                case 'NEU'
                    
                    matA(idIntU,idIntP) = matA(idIntU,idIntP) +                      nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + eta * nx           * nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + eta * ny           * nx * matMel;
                    
                    matA(idIntV,idIntP) = matA(idIntV,idIntP) +                      ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + eta * nx           * ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + eta * ny           * ny * matMel;
                    
                    gnu = nx*rhsUel + ny*rhsVel;
                    rhsA(idIntP) = rhsA(idIntP) - gnu;
                    rhsA(idIntU) = rhsA(idIntU) + eta * gnu * nx;
                    rhsA(idIntV) = rhsA(idIntV) + eta * gnu * ny;
                    
                case 'ABC'
                    
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5                 * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5 * eta * nx      * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + 0.5 * eta * ny      * matMel;
                    
                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5            * nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5 * eta * nx * nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5 * eta * ny * nx * matMel;
                    
                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + 0.5            * ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5 * eta * nx * ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5 * eta * ny * ny * matMel;
                    
                    gchar = rhsPel - etaNeigh * (nx*rhsUel + ny*rhsVel);
                    rhsA(idIntP) = rhsA(idIntP) + gchar * 0.5;
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