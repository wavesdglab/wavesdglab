function [solA, matA, rhsA] = computeSolNum2D_DG1(mesh, dofm, tau)
disp(['--- CALL computeSolNum2D_DG1']);

global k

numDofTRI = dofm.numDofTRI;

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

[matM, ~, matDX, matDY, rhsP] = buildMatrixGlo2D_DG(mesh, dofm);

matA = [
    -1i*k*matM  -matDX                       -matDY                      ;
    -matDX      -1i*k*matM                   sparse(numDofTRI,numDofTRI) ;
    -matDY      sparse(numDofTRI,numDofTRI)  -1i*k*matM                  ];

rhsA = [
    -1/(1i*k)*rhsP ;
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
        
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (triNeigh > 0)
            
            % Get global ID for exterior unknowns
            dofExt = dofLocTriNeigh(facNeigh,:);
            idExtP = dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = idExtP + numDofTRI;
            idExtV = idExtU + numDofTRI;
            
            matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5*tau         * matMel;
            matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5*nx          * matMel;
            matA(idIntP,idIntV) = matA(idIntP,idIntV) + 0.5*ny          * matMel;
            matA(idIntP,idExtP) = matA(idIntP,idExtP) - 0.5*tau         * matMel;
            matA(idIntP,idExtU) = matA(idIntP,idExtU) + 0.5*nx          * matMel;
            matA(idIntP,idExtV) = matA(idIntP,idExtV) + 0.5*ny          * matMel;
            
            matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5*nx          * matMel;
            matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5/tau * nx*nx * matMel;
            matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5/tau * nx*ny * matMel;
            matA(idIntU,idExtP) = matA(idIntU,idExtP) + 0.5*nx          * matMel;
            matA(idIntU,idExtU) = matA(idIntU,idExtU) - 0.5/tau * nx*nx * matMel;
            matA(idIntU,idExtV) = matA(idIntU,idExtV) - 0.5/tau * nx*ny * matMel;
            
            matA(idIntV,idIntP) = matA(idIntV,idIntP) + 0.5*ny          * matMel;
            matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5/tau * nx*ny * matMel;
            matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5/tau * ny*ny * matMel;
            matA(idIntV,idExtP) = matA(idIntV,idExtP) + 0.5*ny          * matMel;
            matA(idIntV,idExtU) = matA(idIntV,idExtU) - 0.5/tau * nx*ny * matMel;
            matA(idIntV,idExtV) = matA(idIntV,idExtV) - 0.5/tau * ny*ny * matMel;
            
        else
            
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            switch tagToBC(mesh.tagEdg(edgGlo))
                case 'DIR'
                    
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + tau         * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + nx          * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + ny          * matMel;
                    
                    gp = rhsPel;
                    rhsA(idIntP) = rhsA(idIntP) + gp * tau;
                    rhsA(idIntU) = rhsA(idIntU) - gp * nx;
                    rhsA(idIntV) = rhsA(idIntV) - gp * ny;
                    
                case 'NEU'
                    
                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + 1      * nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + nx/tau * nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + ny/tau * nx * matMel;
                    
                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + 1      * ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + nx/tau * ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + ny/tau * ny * matMel;
                    
                    gnu = nx*rhsUel + ny*rhsVel;
                    rhsA(idIntP) = rhsA(idIntP) - gnu;
                    rhsA(idIntU) = rhsA(idIntU) + gnu * nx/tau;
                    rhsA(idIntV) = rhsA(idIntV) + gnu * ny/tau;
                    
                case 'ABC'
                    
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) + 0.5           * matMel;
                    matA(idIntP,idIntU) = matA(idIntP,idIntU) + 0.5 * nx      * matMel;
                    matA(idIntP,idIntV) = matA(idIntP,idIntV) + 0.5 * ny      * matMel;
                    
                    matA(idIntU,idIntP) = matA(idIntU,idIntP) + 0.5      * nx * matMel;
                    matA(idIntU,idIntU) = matA(idIntU,idIntU) + 0.5 * nx * nx * matMel;
                    matA(idIntU,idIntV) = matA(idIntU,idIntV) + 0.5 * ny * nx * matMel;
                    
                    matA(idIntV,idIntP) = matA(idIntV,idIntP) + 0.5      * ny * matMel;
                    matA(idIntV,idIntU) = matA(idIntV,idIntU) + 0.5 * nx * ny * matMel;
                    matA(idIntV,idIntV) = matA(idIntV,idIntV) + 0.5 * ny * ny * matMel;
                    
                    gchar = rhsPel - (nx*rhsUel + ny*rhsVel);
                    rhsA(idIntP) = rhsA(idIntP) + gchar * 0.5;
                    rhsA(idIntU) = rhsA(idIntU) - gchar * 0.5 * nx;
                    rhsA(idIntV) = rhsA(idIntV) - gchar * 0.5 * ny;
                    
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