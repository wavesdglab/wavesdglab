function [solA, sysA, condLoc] = computeSolNum2D_UDG(mesh, dofm, tau, prec)

global k

numDofTRI = dofm.numDofTRI;
numDofFAC = dofm.numDofFAC;
numDofPerTRI = dofm.numDofPerTRI;

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);
weightsLinQ = sparse(1:size(weightsLinQ,1), 1:size(weightsLinQ,1), weightsLinQ);
weightsTriQ = sparse(1:size(weightsTriQ,1), 1:size(weightsTriQ,1), weightsTriQ);

% Shape functions
shapeLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeTriDuQ, shapeTriDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% Global matrices
% matII    = sparse(3*numDofTRI, 3*numDofTRI);
% matIIinv = sparse(3*numDofTRI, 3*numDofTRI);
% matIG    = sparse(3*numDofTRI, numDofFAC);
% matGI    = sparse(numDofFAC, 3*numDofTRI);
% matGG    = sparse(numDofFAC, numDofFAC);
% matGGinv = sparse(numDofFAC, numDofFAC);

matIIx = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);
matIIy = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);
matIIv = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);
matIGx = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);
matIGy = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);
matIGv = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);
matGIx = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);
matGIy = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);
matGIv = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);
matGGx = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);
matGGy = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);
matGGv = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);
matIIvInv = zeros(mesh.numTri*3*dofm.numDofPerTRI,3*dofm.numDofPerTRI);
matGGvInv = zeros(mesh.numTri*3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);

rhsI = zeros(3*numDofTRI,1);
rhsG = zeros(numDofFAC,1);

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
    
    % -------------------------------------------------------------------------
    % Surface terms
    % -------------------------------------------------------------------------
    
    matIGel = zeros(3*dofm.numDofPerTRI,3*dofm.numDofPerLIN);
    matGIel = zeros(3*dofm.numDofPerLIN,3*dofm.numDofPerTRI);
    matGGel = zeros(3*dofm.numDofPerLIN,3*dofm.numDofPerLIN);
    rhsGel  = zeros(3*dofm.numDofPerLIN,1);
    
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
        [xQ, yQ] = locToGloLIN(uLinQ, V1, V2);
        Jdxdu = norm(V2-V1) * 0.5;  % [ dx/du ]
        
        % Orientation
        orientation = ones(dofm.numDofPerLIN,1);
        if(n1(fac) > n2(fac))
            orientation(3:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        orientation = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation);
        
        % Shape function
        shapeOrQ = shapeLinQ * orientation;
        
        % Elemental matrices
        matMel = shapeOrQ' * weightsLinQ * shapeOrQ * Jdxdu;
        matIel = sparse(1:size(matMel,1),1:size(matMel,2),1);
        
        % Exterior normal
        nx = normal(fac,1);
        ny = normal(fac,2);
        
        % Local ID for interior unknowns and incoming characteristics
        idLocP = 0*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocU = 1*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocV = 2*dofm.numDofPerTRI + dofm.locFac(fac,:);
        idLocG = (1:dofm.numDofPerLIN) + (fac-1)*dofm.numDofPerLIN;
        
        % Interior contributions
        matIIel(idLocP,idLocP) = matIIel(idLocP,idLocP) + 0.5*tau           * matMel;
        matIIel(idLocP,idLocU) = matIIel(idLocP,idLocU) + 0.5     * nx      * matMel;
        matIIel(idLocP,idLocV) = matIIel(idLocP,idLocV) + 0.5     * ny      * matMel;
        matIGel(idLocP,idLocG) = matIGel(idLocP,idLocG) - 0.5               * matMel;
        
        matIIel(idLocU,idLocP) = matIIel(idLocU,idLocP) + 0.5          * nx * matMel;
        matIIel(idLocU,idLocU) = matIIel(idLocU,idLocU) + 0.5/tau * nx * nx * matMel;
        matIIel(idLocU,idLocV) = matIIel(idLocU,idLocV) + 0.5/tau * nx * ny * matMel;
        matIGel(idLocU,idLocG) = matIGel(idLocU,idLocG) + 0.5/tau      * nx * matMel;
        
        matIIel(idLocV,idLocP) = matIIel(idLocV,idLocP) + 0.5          * ny * matMel;
        matIIel(idLocV,idLocU) = matIIel(idLocV,idLocU) + 0.5/tau * nx * ny * matMel;
        matIIel(idLocV,idLocV) = matIIel(idLocV,idLocV) + 0.5/tau * ny * ny * matMel;
        matIGel(idLocV,idLocG) = matIGel(idLocV,idLocG) + 0.5/tau      * ny * matMel;
        
        % Infos on neighboring element
        triNeigh = mesh.mapTriToTri(tri,fac);
        facNeigh = mesh.mapTriToFac(tri,fac);
        
        if (prec == 0)
            matGGel(idLocG,idLocG) = matMel;
        else
            matGGel(idLocG,idLocG) = matIel;
        end
        
        if (triNeigh > 0)
            
            % Get global ID for exterior unknowns
            idIncG = dofm.locToGloFAC(tri,idLocG);
            dofExt = dofm.locFacNeigh(facNeigh,:);
            idExtP = 0*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtU = 1*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            idExtV = 2*numDofTRI + dofm.locToGloTRI(triNeigh,dofExt);
            
            % matGI(idIncG,idExtP) = -tau * matMel;
            % matGI(idIncG,idExtU) =  nx  * matMel;
            % matGI(idIncG,idExtV) =  ny  * matMel;
            
            idLIN = (tri-1)*3*dofm.numDofPerLIN + (fac-1)*dofm.numDofPerLIN + (1:dofm.numDofPerLIN);
            matGIx(idLIN,:) = idIncG'*ones(1,size([idExtP idExtU idExtV],2));
            matGIy(idLIN,:) = ones(size(idIncG,2),1)*[idExtP idExtU idExtV];
            if (prec == 0)
                matGIv(idLIN,:) = [-tau*matMel nx*matMel ny*matMel];
            else
                matGIv(idLIN,:) = [-tau*matIel nx*matIel ny*matIel];
            end
            
        else
            
            % Solution function
            [solQ, solDxQ, solDyQ, ~] = mySol(xQ, yQ);
            
            % Elemental matrices
            rhsPel = shapeOrQ' * weightsLinQ * solQ * Jdxdu;
            rhsUel = shapeOrQ' * weightsLinQ * solDxQ * Jdxdu / (1i*k);
            rhsVel = shapeOrQ' * weightsLinQ * solDyQ * Jdxdu / (1i*k);
            
            edgGlo = abs(mesh.mapTriToEdg(tri,fac));
            
            if (prec == 0)
                switch tagToBC(mesh.tagEdg(edgGlo))
                    case 'DIR'
                        matGIel(idLocG,idLocP) = + tau * matMel;
                        matGIel(idLocG,idLocU) = + nx  * matMel;
                        matGIel(idLocG,idLocV) = + ny  * matMel;
                        rhsGel(idLocG) = + 2*tau*rhsPel;
                    case 'NEU'
                        matGIel(idLocG,idLocP) = - tau * matMel;
                        matGIel(idLocG,idLocU) = - nx  * matMel;
                        matGIel(idLocG,idLocV) = - ny  * matMel;
                        rhsGel(idLocG) = - 2*(nx*rhsUel + ny*rhsVel);
                    case 'ABC'
                        matGIel(idLocG,idLocP) = + tau * matMel * (1-tau)/(1+tau);
                        matGIel(idLocG,idLocU) = + nx  * matMel * (1-tau)/(1+tau);
                        matGIel(idLocG,idLocV) = + ny  * matMel * (1-tau)/(1+tau);
                        rhsGel(idLocG) = + (rhsPel - (nx*rhsUel + ny*rhsVel)) * (2*tau)/(1+tau);
                    otherwise
                        warning('Error - Bad BC.')
                end
            else
                switch tagToBC(mesh.tagEdg(edgGlo))
                    case 'DIR'
                        matGIel(idLocG,idLocP) = + tau * matIel;
                        matGIel(idLocG,idLocU) = + nx  * matIel;
                        matGIel(idLocG,idLocV) = + ny  * matIel;
                        rhsGel(idLocG) = + 2*tau*matMel\rhsPel;
                    case 'NEU'
                        matGIel(idLocG,idLocP) = - tau * matIel;
                        matGIel(idLocG,idLocU) = - nx  * matIel;
                        matGIel(idLocG,idLocV) = - ny  * matIel;
                        rhsGel(idLocG) = - 2*matMel\(nx*rhsUel + ny*rhsVel);
                    case 'ABC'
                        matGIel(idLocG,idLocP) = + tau * matIel * (1-tau)/(1+tau);
                        matGIel(idLocG,idLocU) = + nx  * matIel * (1-tau)/(1+tau);
                        matGIel(idLocG,idLocV) = + ny  * matIel * (1-tau)/(1+tau);
                        rhsGel(idLocG) = + matMel\(rhsPel - (nx*rhsUel + ny*rhsVel)) * (2*tau)/(1+tau);
                    otherwise
                        warning('Error - Bad BC.')
                end
            end
            
            idLIN = (tri-1)*3*dofm.numDofPerLIN + (fac-1)*dofm.numDofPerLIN + (1:dofm.numDofPerLIN);
            idGloP = 0*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloU = 1*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloV = 2*numDofTRI + dofm.locToGloTRI(tri,idLocP);
            idGloG = dofm.locToGloFAC(tri,idLocG);
            matGIx(idLIN,:) = idGloG'*ones(1,size([idGloP idGloU idGloV],2));
            matGIy(idLIN,:) = ones(size(idGloG,2),1)*[idGloP idGloU idGloV];
            matGIv(idLIN,:) = matGIel(idLocG,[idLocP idLocU idLocV]);
        end
    end
    
    % -------------------------------------------------------------------------
    % Matrix assembling
    % -------------------------------------------------------------------------
    
    dofGloP = 0*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloU = 1*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloV = 2*numDofTRI + dofm.locToGloTRI(tri,:);
    dofGloI = [dofGloP dofGloU dofGloV];
    dofGloG = dofm.locToGloFAC(tri,:);
    
    % matIIinv(dofGloI,dofGloI) = inv(matIIel);
    % matII(dofGloI,dofGloI) = matIIel;
    % matIG(dofGloI,dofGloG) = matIGel;
    % matGI(dofGloG,dofGloI) = matGIel;
    % matGGinv(dofGloG,dofGloG) = inv(matGGel);
    % matGG(dofGloG,dofGloG) = matGGel;
    
    idLIN = (tri-1)*3*dofm.numDofPerLIN + (1:3*dofm.numDofPerLIN);
    idTRI = (tri-1)*3*dofm.numDofPerTRI + (1:3*dofm.numDofPerTRI);
    matIIx(idTRI,:) = dofGloI'*ones(1,size(dofGloI,2));
    matIGx(idTRI,:) = dofGloI'*ones(1,size(dofGloG,2));
    matGGx(idLIN,:) = dofGloG'*ones(1,size(dofGloG,2));
    matIIy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloI;
    matIGy(idTRI,:) = ones(size(dofGloI,2),1)*dofGloG;
    matGGy(idLIN,:) = ones(size(dofGloG,2),1)*dofGloG;
    matIIv(idTRI,:) = matIIel;
    matIGv(idTRI,:) = matIGel;
    matGGv(idLIN,:) = matGGel;
    matIIvInv(idTRI,:) = inv(matIIel);
    matGGvInv(idLIN,:) = inv(matGGel);
    
    condLoc(tri) = cond(full(matIIel));
    
%     figure(2);
%     hold off;
%     diag(matIIel)
%     LLL = size(full(matIIel),1);
%     tmp1 = zeros(LLL,1);
%     tmp2 = zeros(LLL,1);
%     for i=1:LLL
%         tmp1(i) = abs(matIIel(i,i));
%         for j=1:LLL
%             if(i ~= j)
%                 tmp2(i) = tmp2(i) + abs(matIIel(i,j));
%             end
%         end
%     end
%     eigval = eigs(full(matIIel),LLL);
%     min(real(eigval))
%     [tmp1 tmp2]
%     scatter(real(eigval),imag(eigval),'DisplayName','Eigenvalues');
%     pause;
    
    rhsI(dofGloI) = rhsIel;
    rhsG(dofGloG) = rhsGel;
    
end

matII    = sparse(matIIx, matIIy, matIIv, 3*numDofTRI, 3*numDofTRI);
matIG    = sparse(matIGx, matIGy, matIGv, 3*numDofTRI, numDofFAC);
matGI    = sparse(matGIx, matGIy, matGIv, numDofFAC, 3*numDofTRI);
matGG    = sparse(matGGx, matGGy, matGGv, numDofFAC, numDofFAC);
matIIinv = sparse(matIIx, matIIy, matIIvInv, 3*numDofTRI, 3*numDofTRI);
matGGinv = sparse(matGGx, matGGy, matGGvInv, numDofFAC, numDofFAC);

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

matS = matGG - matGI*(matIIinv*matIG);
rhsS = rhsG - matGI*(matIIinv*rhsI);
solG = matS\rhsS;
solI = matIIinv*(rhsI-matIG*solG);
solA = [ solI ; solG ];

matPhy = matII - matIG*(matGGinv*matGI);
rhsPhy = rhsI - matIG*(matGGinv*rhsG);

sysA.matII = matII;
sysA.matIG = matIG;
sysA.matGI = matGI;
sysA.matGG = matGG;
sysA.matIIinv = matIIinv;
sysA.matGGinv = matGGinv;

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
        warning('Error - No valid BC has been set on the East.')
end
end