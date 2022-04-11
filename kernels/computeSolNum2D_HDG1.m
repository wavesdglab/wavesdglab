function [solA, matA, rhsA] = computeSolNum2D_HDG1(mesh, dofm, tau)

fprintf('Solver  : Call computeSolNumHDG1\n');

global k

numDofTRI = dofm.numDofTRI;
numDofLIN = dofm.numDofLIN;

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

% Quadrature
degreeQ = 4*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);

% Shape functions (f, dfdu, dfdv)
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);
[shapeDuQ, shapeDvQ] = functionsShapeDerTRI(uQ, vQ, dofm.degree);

% Global matrices
matM = sparse(dofm.numDofTRI, dofm.numDofTRI);
matDX = sparse(dofm.numDofTRI, dofm.numDofTRI);
matDY = sparse(dofm.numDofTRI, dofm.numDofTRI);
rhsP = zeros(dofm.numDofTRI, 1);

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
    
    % RHS function
    [~, ~, ~, rhsQ] = mySol(xQ, yQ);
    
    % Orientation
    orientation = ones(dofm.numDofPerTRI,1);
    if(ver(1) > ver(2))
        orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(ver(2) > ver(3))
        orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(ver(3) > ver(1))
        orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    orientation = sparse(1:dofm.numDofPerTRI, 1:dofm.numDofPerTRI, orientation);
    
    % Shape functions (f, dfdx, dfdy) with orientation
    shapeOrQ = shapeQ * orientation;
    shapeDxQ = (shapeDuQ * Jdudx(1,1) + shapeDvQ * Jdudx(2,1)) * orientation;
    shapeDyQ = (shapeDuQ * Jdudx(1,2) + shapeDvQ * Jdudx(2,2)) * orientation;
    
    % Elemental matrices
    matMel = shapeOrQ' * weights * shapeOrQ * detJdxdu;
    matDXel = shapeDxQ' * weights * shapeOrQ * detJdxdu;
    matDYel = shapeDyQ' * weights * shapeOrQ * detJdxdu;
    rhsPel = shapeOrQ' * weights * rhsQ * detJdxdu;
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matM(dof,dof) = matMel;
    matDX(dof,dof) = matDXel;
    matDY(dof,dof) = matDYel;
    rhsP(dof) = rhsPel;
    
end

matA = [
    -1i*k*matM                   -matDX                       -matDY                       sparse(numDofTRI,numDofLIN) ;
    -matDX                       -1i*k*matM                   sparse(numDofTRI,numDofTRI)  sparse(numDofTRI,numDofLIN) ;
    -matDY                       sparse(numDofTRI,numDofTRI)  -1i*k*matM                   sparse(numDofTRI,numDofLIN) ;
    sparse(numDofLIN,numDofTRI)  sparse(numDofLIN,numDofTRI)  sparse(numDofLIN,numDofTRI)  sparse(numDofLIN,numDofLIN) ];

rhsA = [
    -1/(1i*k)*rhsP ;
    zeros(numDofTRI,1)  ;
    zeros(numDofTRI,1)  ;
    zeros(numDofLIN,1)  ];

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
        
        % Global ID for edge unknowns
        edgGlo = abs(mesh.mapTriToEdg(tri,fac));
        idIntS = 3*numDofTRI + dofm.locToGloLIN(edgGlo,:);
        if(mesh.mapTriToEdg(tri,fac) < 0)
            tmp = idIntS;
            idIntS(1) = tmp(2);
            idIntS(2) = tmp(1);
        end
        
        % Mapping
        V1 = mesh.coord(n1(fac),:);
        V2 = mesh.coord(n2(fac),:);
        [xQ, yQ] = locToGloLIN(uQ, V1, V2);
        Jdxdu = norm(V2-V1) * 0.5;  % [ dx/du ]
        
        % Solution function
        [solQ, solDxQ, solDyQ, ~] = mySol(xQ, yQ);
        
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
        rhsUel = shapeOrQ' * weights * solDxQ * Jdxdu / (1i*k);
        rhsVel = shapeOrQ' * weights * solDyQ * Jdxdu / (1i*k);
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        % Surface terms for the volume fields
        matA(idIntP,idIntP) = matA(idIntP,idIntP) + tau * matMel;
        matA(idIntP,idIntU) = matA(idIntP,idIntU) + nx  * matMel;
        matA(idIntP,idIntV) = matA(idIntP,idIntV) + ny  * matMel;
        matA(idIntP,idIntS) = matA(idIntP,idIntS) - tau * matMel;
        matA(idIntU,idIntS) = matA(idIntU,idIntS) + nx  * matMel;
        matA(idIntV,idIntS) = matA(idIntV,idIntS) + ny  * matMel;
        
        % Surface terms for the surface field
        triNeigh = mesh.mapTriToTri(tri,fac);
        if (triNeigh > 0)
            matA(idIntS,idIntP) = matA(idIntS,idIntP) + tau * matMel;
            matA(idIntS,idIntU) = matA(idIntS,idIntU) + nx  * matMel;
            matA(idIntS,idIntV) = matA(idIntS,idIntV) + ny  * matMel;
            matA(idIntS,idIntS) = matA(idIntS,idIntS) - tau * matMel;
        else
            
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            switch tagToBC(mesh.tagEdg(edgGlo))
                case 'DIR'
                    matA(idIntS,idIntS) = matA(idIntS,idIntS) + matMel;
                    rhsA(idIntS) = rhsA(idIntS) + rhsPel;
                case 'NEU'
                    matA(idIntS,idIntP) = matA(idIntS,idIntP) + tau * matMel;
                    matA(idIntS,idIntU) = matA(idIntS,idIntU) + nx  * matMel;
                    matA(idIntS,idIntV) = matA(idIntS,idIntV) + ny  * matMel;
                    matA(idIntS,idIntS) = matA(idIntS,idIntS) - tau * matMel;
                    rhsA(idIntS) = rhsA(idIntS) + (nx*rhsUel + ny*rhsVel);
                case 'ABC'
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + (1-tau) * matMel;
                    matA(idIntP,idIntS) = matA(idIntP,idIntS) - (1-tau) * matMel;
                    % 
                    matA(idIntS,idIntP) = matA(idIntS,idIntP) +       matMel;
                    matA(idIntS,idIntU) = matA(idIntS,idIntU) + nx  * matMel;
                    matA(idIntS,idIntV) = matA(idIntS,idIntV) + ny  * matMel;
                    matA(idIntS,idIntS) = matA(idIntS,idIntS) - 2   * matMel;
                    rhsA(idIntS) = rhsA(idIntS) + (nx*rhsUel + ny*rhsVel - rhsPel);
                otherwise
                    warning('Error - Bad BC.');
            end
        end
    end
end

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

solA = matA\rhsA;
solA = solA(1:numDofTRI);

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
        warning('Error - No valid BC has been set on the East.')
end
end