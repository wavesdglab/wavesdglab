function [solA, matA, rhsA] = computeSolNum2D_UDG1(mesh, dofm, tau)
disp(['--- CALL computeSolNum2D_UDG1']);

global k

numDofTRI = dofm.numDofTRI;
numDofFAC = dofm.numDofFAC;

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

[matM, ~, matDX, matDY, rhsP] = buildMatrixGlo2D_DG(mesh, dofm);

matA = [
    -1i*k*matM                   -matDX                       -matDY                       sparse(numDofTRI,numDofFAC) ;
    -matDX                       -1i*k*matM                   sparse(numDofTRI,numDofTRI)  sparse(numDofTRI,numDofFAC) ;
    -matDY                       sparse(numDofTRI,numDofTRI)  -1i*k*matM                   sparse(numDofTRI,numDofFAC) ;
    sparse(numDofFAC,numDofTRI)  sparse(numDofFAC,numDofTRI)  sparse(numDofFAC,numDofTRI)  sparse(numDofFAC,numDofFAC) ];

rhsA = [
    -1/(1i*k)*rhsP ;
    zeros(numDofTRI,1) ;
    zeros(numDofTRI,1) ;
    zeros(numDofFAC,1) ];

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

% Quadrature
degreeQ = 2*dofm.degree;
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
        
        % Global ID for interior unknowns and incoming characteristics
        dofInt = dofm.locFac(fac,:);
        idIntP = 0*numDofTRI + dofm.locToGloTRI(tri,dofInt);
        idIntU = 1*numDofTRI + dofm.locToGloTRI(tri,dofInt);
        idIntV = 2*numDofTRI + dofm.locToGloTRI(tri,dofInt);
        idIncG = 3*numDofTRI + dofm.locToGloFAC(tri,(1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN);
        
        % Interior contributions
        matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau           * matMel;
        matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5     * nx      * matMel;
        matA(idIntP,idIntV) = matA(idIntP,idIntV) + 0.5     * ny      * matMel;
        matA(idIntP,idIncG) = matA(idIntP,idIncG) - 0.5               * matMel;
        
        matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5          * nx * matMel;
        matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau * nx * nx * matMel;
        matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5/tau * nx * ny * matMel;
        matA(idIntU,idIncG) = matA(idIntU,idIncG) + 0.5/tau      * nx * matMel;
        
        matA(idIntV,idIntP) = matA(idIntV,idIntP) + 0.5          * ny * matMel;
        matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5/tau * nx * ny * matMel;
        matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5/tau * ny * ny * matMel;
        matA(idIntV,idIncG) = matA(idIntV,idIncG) + 0.5/tau      * ny * matMel;
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        matA(idIncG,idIncG) = matMel;
        
        if (triNeigh > 0)
            
            % Get global ID for exterior unknowns
            dofExt = dofm.locFacNeigh(facNeigh,:);
            idExtP = 0*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = 1*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtV = 2*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            
            matA(idIncG,idExtP) = matA(idIncG,idExtP) - tau * matMel;
            matA(idIncG,idExtU) = matA(idIncG,idExtU) + nx  * matMel;
            matA(idIncG,idExtV) = matA(idIncG,idExtV) + ny  * matMel;
            
        else
            
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            switch tagToBC(mesh.tagEdg(edgGlo))
                case 'DIR'
                    matA(idIncG,idIntP) = matA(idIncG,idIntP) + tau * matMel;
                    matA(idIncG,idIntU) = matA(idIncG,idIntU) + nx  * matMel;
                    matA(idIncG,idIntV) = matA(idIncG,idIntV) + ny  * matMel;
                    rhsA(idIncG) = rhsA(idIncG) + 2*tau*rhsPel;
                case 'NEU'
                    matA(idIncG,idIntP) = matA(idIncG,idIntP) - tau * matMel;
                    matA(idIncG,idIntU) = matA(idIncG,idIntU) - nx  * matMel;
                    matA(idIncG,idIntV) = matA(idIncG,idIntV) - ny  * matMel;
                    rhsA(idIncG) = rhsA(idIncG) - 2*(nx*rhsUel + ny*rhsVel);
                case 'ABC'
                    
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*(1-tau)           * matMel;
                    
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5*(1-1/tau) * nx * nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5*(1-1/tau) * nx * ny * matMel;
                    matA(idIntU,idIncG) = matA(idIntU,idIncG) + 0.5*(1-1/tau)      * nx * matMel;
                    
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5*(1-1/tau) * nx * ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5*(1-1/tau) * ny * ny * matMel;
                    matA(idIntV,idIncG) = matA(idIntV,idIncG) + 0.5*(1-1/tau)      * ny * matMel;
                    % 
                    rhsA(idIncG) = rhsA(idIncG) + (tau*rhsPel - (nx*rhsUel + ny*rhsVel));
                otherwise
                    warning('Error - Bad BC.')
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