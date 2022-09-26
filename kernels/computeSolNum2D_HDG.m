function [solA, sysA, condLoc] = computeSolNum2D_HDG(mesh, dofm, tau, prec)

global k

numDofTRI = dofm.numDofTRI;
numDofLIN = dofm.numDofLIN;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);
weightsLinQ = sparse(1:size(weightsLinQ,1), 1:size(weightsLinQ,1), weightsLinQ);
weightsTriQ = sparse(1:size(weightsTriQ,1), 1:size(weightsTriQ,1), weightsTriQ);

% Shape functions
shapeIntQ = functionsShapeINT(uLinQ, dofm.degree);
shapeLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeTriDuQ, shapeTriDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% Global matrices
% matII    = sparse(3*numDofTRI, 3*numDofTRI);
% matIG    = sparse(3*numDofTRI, numDofLIN);
% matGI    = sparse(numDofLIN, 3*numDofTRI);
% matGG    = sparse(numDofLIN, numDofLIN);
% matIIinv = sparse(3*numDofTRI, 3*numDofTRI);
rhsI     = zeros(3*numDofTRI,1);
rhsG     = zeros(numDofLIN,1);

matIIx = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);
matIIy = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);
matIIv = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);
matIGx = zeros(3*dofm.numDofPerTRI,mesh.numTri*3*dofm.numDofPerLIN);
matIGy = zeros(3*dofm.numDofPerTRI,mesh.numTri*3*dofm.numDofPerLIN);
matIGv = zeros(3*dofm.numDofPerTRI,mesh.numTri*3*dofm.numDofPerLIN);
matGIx = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerTRI);
matGIy = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerTRI);
matGIv = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerTRI);
matGGx = zeros(numDofLIN,dofm.numDofPerLIN);
matGGy = zeros(numDofLIN,dofm.numDofPerLIN);
matGGv = zeros(numDofLIN,dofm.numDofPerLIN);
matIIvInv = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);

condLoc = zeros(mesh.numTri,1);

for tri=1:mesh.numTri
    
    % -------------------------------------------------------------------------
    % Volume terms
    % -------------------------------------------------------------------------
    
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
    
    % Shape functions (f, dfdx, dfdy) with orientation
    shapeOrQ = shapeTriQ * orientation;
    shapeDxQ = (shapeTriDuQ * Jdudx(1,1) + shapeTriDvQ * Jdudx(2,1)) * orientation;
    shapeDyQ = (shapeTriDuQ * Jdudx(1,2) + shapeTriDvQ * Jdudx(2,2)) * orientation;
    
    % RHS function
    [~, ~, ~, rhsQ] = mySol(xQ, yQ);
    
    % Elemental matrices
    matMel = shapeOrQ' * weightsTriQ * shapeOrQ * detJdxdu;
    matDXel = shapeDxQ' * weightsTriQ * shapeOrQ * detJdxdu;
    matDYel = shapeDyQ' * weightsTriQ * shapeOrQ * detJdxdu;
    vecRHSel = shapeOrQ' * weightsTriQ * rhsQ * detJdxdu;
    
    matIIel = [
        -1i*k*matMel  -matDXel                           -matDYel                          ;
        -matDXel      -1i*k*matMel                       sparse(numDofPerTRI,numDofPerTRI) ;
        -matDYel      sparse(numDofPerTRI,numDofPerTRI)  -1i*k*matMel                      ];
    
    rhsIel = [
        -1/(1i*k)*vecRHSel    ;
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
        
        matIGel = zeros(3*dofm.numDofPerTRI,dofm.numDofPerLIN);
        matGIel = zeros(dofm.numDofPerLIN,3*dofm.numDofPerTRI);
        matGGel = zeros(dofm.numDofPerLIN,dofm.numDofPerLIN);
        rhsGel = zeros(dofm.numDofPerLIN,1);
        
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
        
        orientation2 = ones(dofm.numDofPerLIN,1);
        if(n1(fac) > n2(fac))
            orientation2(1:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerLIN-1);
        end
        orientation2 = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation2);
        
        % Shape function
        shapeInQ = shapeIntQ * orientation2;
        shapeOrQ = shapeLinQ * orientation;
        
        % Elemental matrices
        matMel = shapeOrQ' * weightsLinQ * shapeOrQ * Jdxdu;
        matLel = shapeOrQ' * weightsLinQ * shapeInQ * Jdxdu;
        matIel = sparse(1:size(matMel,1),1:size(matMel,2),1);
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        % Local ID for interior unknowns and incoming characteristics
        idLocP = 0*numDofPerTRI + dofm.locFac(fac,:);
        idLocU = 1*numDofPerTRI + dofm.locFac(fac,:);
        idLocV = 2*numDofPerTRI + dofm.locFac(fac,:);
        
        % Surface terms for the volume fields
        matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + tau * matMel;
        matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + nx  * matMel;
        matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + ny  * matMel;
        if (prec == 10)
            matIGel(idLocP,:) = -tau * matLel;
            matIGel(idLocU,:) =  nx  * matLel;
            matIGel(idLocV,:) =  ny  * matLel;
        else
            matIGel(idLocP,:) = -tau * matMel;
            matIGel(idLocU,:) =  nx  * matMel;
            matIGel(idLocV,:) =  ny  * matMel;
        end
        
        % Surface terms for the surface field
        triNeigh = mesh.mapTriToTri(tri,fac);
        if (triNeigh > 0)
            
            if (prec == 0)
                matGIel(:,idLocP) = tau * matMel;
                matGIel(:,idLocU) = nx  * matMel;
                matGIel(:,idLocV) = ny  * matMel;
                matGGel = -tau * matMel;
            elseif (prec == 1)
                matGIel(:,idLocP) = tau * matIel;
                matGIel(:,idLocU) = nx  * matIel;
                matGIel(:,idLocV) = ny  * matIel;
                matGGel = -tau * matIel;
            elseif (prec == 2)
                matGIel(:,idLocP) = tau * matMel;
                matGIel(:,idLocU) = nx  * matMel;
                matGIel(:,idLocV) = ny  * matMel;
                matGGel = -tau * matMel;
            elseif (prec == 10)
                matGIel(:,idLocP) = -0.5 * tau * matLel'/tau/Jdxdu;
                matGIel(:,idLocU) = -0.5 * nx  * matLel'/tau/Jdxdu;
                matGIel(:,idLocV) = -0.5 * ny  * matLel'/tau/Jdxdu;
                matGGel = 0.5 * matIel;
            end
            
        else
            
            % Solution function
            [solQ, solDxQ, solDyQ, ~] = mySol(xQ, yQ);
            rhsPel = shapeOrQ' * weightsLinQ * solQ * Jdxdu;
            rhsUel = shapeOrQ' * weightsLinQ * solDxQ * Jdxdu / (1i*k);
            rhsVel = shapeOrQ' * weightsLinQ * solDyQ * Jdxdu / (1i*k);
            
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            if (prec == 0)
                switch tagToBC(mesh.tagEdg(edgGlo))
                    case 'DIR'
                        matGGel = matMel;
                        rhsGel = rhsPel;
                    case 'NEU'
                        matGIel(:,idLocP) = tau * matMel;
                        matGIel(:,idLocU) = nx  * matMel;
                        matGIel(:,idLocV) = ny  * matMel;
                        matGGel = -tau * matMel;
                        rhsGel = nx*rhsUel + ny*rhsVel;
                    case 'ABC'
                        matGIel(:,idLocP) = -tau * matMel;
                        matGIel(:,idLocU) = -nx  * matMel;
                        matGIel(:,idLocV) = -ny  * matMel;
                        matGGel = (1+tau) * matMel;
                        rhsGel  = rhsPel - nx*rhsUel - ny*rhsVel;
                    otherwise
                        warning('Error - Bad BC.');
                end
            elseif (prec == 1)
                switch tagToBC(mesh.tagEdg(edgGlo))
                    case 'DIR'
                        matGGel = matIel;
                        rhsGel = matMel\rhsPel;
                    case 'NEU'
                        matGIel(:,idLocP) = tau * matIel;
                        matGIel(:,idLocU) = nx  * matIel;
                        matGIel(:,idLocV) = ny  * matIel;
                        matGGel = -tau * matIel;
                        rhsGel = matMel\(nx*rhsUel + ny*rhsVel);
                    case 'ABC'
                        matGIel(:,idLocP) = -tau * matIel;
                        matGIel(:,idLocU) = -nx  * matIel;
                        matGIel(:,idLocV) = -ny  * matIel;
                        matGGel = (1+tau) * matIel;
                        rhsGel  = matMel\(rhsPel - nx*rhsUel - ny*rhsVel);
                    otherwise
                        warning('Error - Bad BC.');
                end
            elseif (prec == 2)
                switch tagToBC(mesh.tagEdg(edgGlo))
                    case 'DIR'
                        matGGel = matMel;
                        rhsGel = rhsPel;
                    case 'NEU'
                        matGIel(:,idLocP) = tau * matMel;
                        matGIel(:,idLocU) = nx  * matMel;
                        matGIel(:,idLocV) = ny  * matMel;
                        matGGel = -tau * matMel;
                        rhsGel = nx*rhsUel + ny*rhsVel;
                    case 'ABC'
                        matGIel(:,idLocP) = -tau * matMel;
                        matGIel(:,idLocU) = -nx  * matMel;
                        matGIel(:,idLocV) = -ny  * matMel;
                        matGGel = (1+tau) * matMel;
                        rhsGel  = rhsPel - nx*rhsUel - ny*rhsVel;
                    otherwise
                        warning('Error - Bad BC.');
                end
            elseif (prec == 10)
                rhsPel = shapeInQ' * weightsLinQ * solQ;
                rhsUel = shapeInQ' * weightsLinQ * solDxQ / (1i*k);
                rhsVel = shapeInQ' * weightsLinQ * solDyQ / (1i*k);
                switch tagToBC(mesh.tagEdg(edgGlo))
                    case 'DIR'
                        matGGel = matIel;
                        rhsGel = rhsPel;
                    case 'NEU'
                        matGIel(:,idLocP) = -tau * matLel' / tau / Jdxdu;
                        matGIel(:,idLocU) = -nx  * matLel' / tau / Jdxdu;
                        matGIel(:,idLocV) = -ny  * matLel' / tau / Jdxdu;
                        matGGel = matIel;
                        rhsGel = -(nx*rhsUel + ny*rhsVel) / tau;
                    case 'ABC'
                        matGIel(:,idLocP) = -tau * matLel' / (1+tau) / Jdxdu;
                        matGIel(:,idLocU) = -nx  * matLel' / (1+tau) / Jdxdu;
                        matGIel(:,idLocV) = -ny  * matLel' / (1+tau) / Jdxdu;
                        matGGel = matIel;
                        rhsGel  = (rhsPel - nx*rhsUel - ny*rhsVel) / (1+tau);
                    otherwise
                        warning('Error - Bad BC.');
                end
            end
        end
        
        % Global ID for edge unknowns
        edgGlo = abs(mesh.mapTriToEdg(tri,fac));
        dofGloG = dofm.locToGloLIN(edgGlo,:);
        dofLocG = 1:dofm.numDofPerLIN;
        if (prec ~= 10)
            if(mesh.mapTriToEdg(tri,fac) < 0)
                dofLocG(1) = 2;
                dofLocG(2) = 1;
            end
        end
        
        % -------------------------------------------------------------------------
        % Matrix assembling
        % -------------------------------------------------------------------------
        
        idLIN = (tri-1)*3*dofm.numDofPerLIN + (fac-1)*dofm.numDofPerLIN + (1:dofm.numDofPerLIN);
        matGIx(idLIN,:) = dofGloG'*ones(1,size(dofGloI,2));
        matGIy(idLIN,:) = ones(size(dofGloG,2),1)*dofGloI;
        matGIv(idLIN,:) = matGIel(dofLocG,:);
        matIGx(:,idLIN) = dofGloI'*ones(1,size(dofGloG,2));
        matIGy(:,idLIN) = ones(size(dofGloI,2),1)*dofGloG;
        matIGv(:,idLIN) = matIGel(:,dofLocG);
        matGGx(dofGloG,:) = dofGloG'*ones(1,size(dofGloG,2));
        matGGy(dofGloG,:) = ones(size(dofGloG,2),1)*dofGloG;
        matGGv(dofGloG,:) = matGGv(dofGloG,:) + matGGel(dofLocG,dofLocG);
        
        % matIG(dofGloI,dofGloG) = matIGel(:,dofLocG);
        % matGI(dofGloG,dofGloI) = matGIel(dofLocG,:);
        % matGG(dofGloG,dofGloG) = matGG(dofGloG,dofGloG) + matGGel(dofLocG,dofLocG);
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
    
    condLoc(tri) = cond(full(matIIel));
    
    % matIIinv(dofGloI,dofGloI) = inv(matIIel);
    % matII(dofGloI,dofGloI) = matIIel;
    rhsI(dofGloI) = rhsIel;
    
end

matII    = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG    = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofLIN);
matGI    = sparse(matGIx, matGIy, matGIv, numDofLIN, 3*numDofTRI);
matGG    = sparse(matGGx, matGGy, matGGv, numDofLIN, numDofLIN);
matIIinv = sparse(matIIx, matIIy, matIIvInv, 3*numDofTRI, 3*numDofTRI);

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

matS = matGG - matGI*(matIIinv*matIG);
rhsS = rhsG - matGI*(matIIinv*rhsI);

sysA.precL = sparse(1:size(matGG,1),1:size(matGG,1),1);
sysA.precR = sparse(1:size(matGG,1),1:size(matGG,1),1);
if (prec == 2)
    [eigenvecGG,eigenvalGG] = eigs(matGG,size(matGG,1));
    sysA.precL = sqrt(eigenvalGG)\eigenvecGG';
    sysA.precR = eigenvecGG/sqrt(eigenvalGG);
end
matS = sysA.precL*matS*sysA.precR;
rhsS = sysA.precL*rhsS;
solG = matS\rhsS;
solG = sysA.precR*solG;

solI = matIIinv*(rhsI-matIG*solG);
solA = [ solI ; solG ];

matPhy = matII - matIG*(matGG\matGI);
rhsPhy = rhsI - matIG*(matGG\rhsG);

sysA.matII = matII;
sysA.matIG = matIG;
sysA.matGI = matGI;
sysA.matGG = matGG;
sysA.matIIinv = matIIinv;


sysA.rhsI = rhsI;
sysA.rhsG = rhsG;

sysA.matA = [ matII matIG ; matGI matGG ];
sysA.rhsA = [ rhsI ; rhsG ];
sysA.matS = matS;
sysA.rhsS = rhsS;
sysA.matPhy = matPhy;
sysA.rhsPhy = rhsPhy;

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
        warning('Error - No valid BC has been set.')
end
end