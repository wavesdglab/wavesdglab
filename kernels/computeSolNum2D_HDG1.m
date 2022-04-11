function [solA, matA, rhsA] = computeSolNum2D_HDG1(mesh, dofm, tau)
disp(['--- CALL computeSolNum2D_HDG1']);

global k

numDofTRI = dofm.numDofTRI;
numDofLIN = dofm.numDofLIN;

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

[matM, ~, matDX, matDY, rhsP] = buildMatrixGlo2D_DG(mesh, dofm);

matA = [
    -1i*k*matM                   -matDX                       -matDY                       sparse(numDofTRI,numDofLIN) ;
    -matDX                       -1i*k*matM                   sparse(numDofTRI,numDofTRI)  sparse(numDofTRI,numDofLIN) ;
    -matDY                       sparse(numDofTRI,numDofTRI)  -1i*k*matM                   sparse(numDofTRI,numDofLIN) ;
    sparse(numDofLIN,numDofTRI)  sparse(numDofLIN,numDofTRI)  sparse(numDofLIN,numDofTRI)  sparse(numDofLIN,numDofLIN) ];

rhsA = [
    -1/(1i*k)*rhsP ;
    zeros(numDofTRI,1) ;
    zeros(numDofTRI,1) ;
    zeros(numDofLIN,1) ];

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

% Quadrature
degreeQ = 4*dofm.degree;
[uQ, weights] = quadratureGaussLIN(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);

% Shape functions
shapeQ = functionsShapeLIN(uQ, dofm.degree);

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
        
        % Elemental matrices
        matMel = shapeOrQ' * weights * shapeOrQ * Jdxdu;
        rhsPel = shapeOrQ' * weights * solQ * Jdxdu;
        rhsUel = shapeOrQ' * weights * solDxQ * Jdxdu / (1i*k);
        rhsVel = shapeOrQ' * weights * solDyQ * Jdxdu / (1i*k);
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        % Global ID for interior unknowns
        dofInt = dofm.locFac(fac,:);
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