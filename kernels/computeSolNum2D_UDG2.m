function [solA, matA, rhsA] = computeSolNum2D_UDG2(mesh, dofm, tau)
disp(['--- CALL computeSolNum2D_UDG2']);

global k

numDofTRI = dofm.numDofTRI;
numDofFAC = dofm.numDofFAC;

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

[matM, matK, ~, ~, rhsP] = buildMatrixGlo2D_DG(mesh, dofm);

matA = [
    matK - k^2*matM              sparse(numDofTRI,numDofFAC) ;
    sparse(numDofFAC,numDofTRI)  sparse(1:numDofFAC,1:numDofFAC,1) ];

rhsA = [
    rhsP ;
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

dofDIR = [];

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
        rhsDXel = shapeOrQ' * weights * solDxQ * Jdxdu;
        rhsDYel = shapeOrQ' * weights * solDyQ * Jdxdu;
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        % Global ID for interior unknowns and incoming characteristics
        dofInt = dofm.locFac(fac,:);
        idIntP = 0*numDofTRI + dofm.locToGloTRI(tri,dofInt);
        idIncG = 1*numDofTRI + dofm.locToGloFAC(tri,(1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN);
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        matA(idIncG,idIncG) = matMel;
        
        if (triNeigh > 0)
            
            matA(idIntP,idIntP) = matA(idIntP,idIntP) - 1i*tau*k * matMel;
            matA(idIntP,idIncG) = matA(idIntP,idIncG) -            matMel;
            
            % Get global ID for exterior unknowns
            dofExt = dofm.locFacNeigh(facNeigh,:);
            idExtP = 0*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idOutG = 1*numDofTRI + dofm.locToGloFAC(triNeigh,(1:dofm.numDofPerLIN) + (facNeigh-1)*dofm.numDofPerLIN);
            
            matA(idIncG,idOutG) = matA(idIncG,idOutG) +            matMel;
            matA(idIncG,idIntP) = matA(idIncG,idIntP) + 2i*tau*k * matMel;
            
        else
            
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            switch tagToBC(mesh.tagEdg(edgGlo))
                case 'DIR'
                    dofDIR = [dofDIR ; idIntP];
                case 'NEU'
                    rhsA(idIntP) = rhsA(idIntP) + (nx*rhsDXel(idIntP) + ny*rhsDYel(idIntP));
                case 'ABC'
                    matA(idIntP,idIntP) = matA(idIntP,idIntP) - 1i*k * matMel;
                    rhsA(idIntP) = rhsA(idIntP) + (nx*rhsDXel(idIntP) + ny*rhsDYel(idIntP) - 1i*k*rhsPel(idIntP));
                otherwise
                    warning('Error - Bad BC.')
            end
        end
    end
end

if(~isempty(dofDIR))
    solP = computeSolProjL2_2D_DG(mesh, dofm);
    dofDIR = unique(dofDIR);
    rhsA = rhsA - matA(:,dofDIR)*solP(dofDIR);
    rhsA(dofDIR) = solP(dofDIR);
    matA(dofDIR,:) = 0;
    matA(:,dofDIR) = 0;
    matA(dofDIR,dofDIR) = eye(size(dofDIR,1),size(dofDIR,1));
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