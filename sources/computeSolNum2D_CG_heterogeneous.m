% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Authors: Axel Modave, Timothée Raynaud

function [solA, sysA] = computeSolNum2D_CG_heterogeneous2(mesh, dofm, PREC)

global kArray rhoArray
global edgTagToBC
global LdomX LdomY LpmlX LpmlY Rdom Rpml
global pntSouTag pntSouVal

% -------------------------------------------------------------------------
% Quadrature
% -------------------------------------------------------------------------

% Quadrature
degreeQ = 2*dofm.degree;
[uLinQ, weightsLinQ] = quadratureGaussLIN(degreeQ);
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);

% Shape functions
shapeLinQ = functionsShapeLIN(uLinQ, dofm.degree);
shapeTriQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);
[shapeDuQ, shapeDvQ] = functionsShapeDerTRI(uTriQ, vTriQ, dofm.degree);

% -------------------------------------------------------------------------
% Volume terms
% -------------------------------------------------------------------------

matIv  = zeros(mesh.numTri * dofm.numDofPerTRI, dofm.numDofPerTRI);
matJv  = zeros(mesh.numTri * dofm.numDofPerTRI, dofm.numDofPerTRI);
matAv  = zeros(mesh.numTri * dofm.numDofPerTRI, dofm.numDofPerTRI);
% matMv  = zeros(mesh.numTri * dofm.numDofPerTRI, dofm.numDofPerTRI);
% matPv  = zeros(mesh.numTri * dofm.numDofPerTRI, dofm.numDofPerTRI);
rhsA = zeros(dofm.numDofTRI, 1);

tic
myWaitbar = waitbar(0,'   Volume terms...');
for tri=1:mesh.numTri
    waitbar(tri/mesh.numTri,myWaitbar,['   Volume terms... (' int2str(100*tri/mesh.numTri) '%)']);
    
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
    shapeOrQ = shapeTriQ * orientation;
    shapeDxQ = (shapeDuQ * Jdudx(1,1) + shapeDvQ * Jdudx(2,1)) * orientation;
    shapeDyQ = (shapeDuQ * Jdudx(1,2) + shapeDvQ * Jdudx(2,2)) * orientation;
    
    % Coefficients on the element
    k = kArray(tri);
    rho = rhoArray(tri);
    
    % RHS term
    rhsQ = mySourceVolume(xQ, yQ);
    weightsQ = weightsTriQ .* detJdxdu;
    rhsPel = transpose(shapeOrQ) * (weightsQ .* rhsQ);
    
    % Elemental matrices/vectors
    weightsQ = weightsTriQ .* detJdxdu;
    matMel = transpose(shapeOrQ) * (weightsQ .* shapeOrQ);
    matKel = transpose(shapeDxQ) * (weightsQ .* shapeDxQ) + transpose(shapeDyQ) * (weightsQ .* shapeDyQ);
    matAel = (1/rho) * matKel - (k^2/rho) * matMel;
    
    % PML stretching (rectangular PML)
    if(~isempty(LdomX) && ~isempty(LdomY))
        if ((mean(abs(xQ)) >= LdomX) || (mean(abs(yQ)) >= LdomY))
            sigmaPmlX = (LdomX <= abs(xQ))./(LdomX+LpmlX-abs(xQ));
            sigmaPmlY = (LdomY <= abs(yQ))./(LdomY+LpmlY-abs(yQ));
            gammaPmlX = ones(size(xQ)) - sigmaPmlX/(1i*k);
            gammaPmlY = ones(size(yQ)) - sigmaPmlY/(1i*k);
            detJdxdu = detJdxdu .* gammaPmlX .* gammaPmlY;
            shapeDxQ = shapeDxQ ./ gammaPmlX;
            shapeDyQ = shapeDyQ ./ gammaPmlY;
            
            % Elemental matrices in PML
            weightsQ = weightsTriQ .* detJdxdu;
            matMpml = transpose(shapeOrQ) * (weightsQ .* shapeOrQ);
            matKpml = transpose(shapeDxQ) * (weightsQ .* shapeDxQ) + transpose(shapeDyQ) * (weightsQ .* shapeDyQ);
            matApml = (1/rho) * matKpml - (k^2/rho) * matMpml;
            
            % Elemental RHS vector in PML
            rhsProj = matMel\rhsPel;
            rhsPel = matApml * rhsProj - matAel * rhsProj;
            matAel = matApml;
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
            invJacXX = (1./gammaPmlR) .* cosT.*cosT + (1./gammaPmlT) .* (sinT.*sinT);
            invJacXY = (1./gammaPmlR) .* cosT.*sinT - (1./gammaPmlT) .* (cosT.*sinT);
            invJacYY = (1./gammaPmlR) .* sinT.*sinT + (1./gammaPmlT) .* (cosT.*cosT);
            detJdxdu = detJdxdu .* gammaPmlR .* gammaPmlT;
            shapeDxQnew = invJacXX .* shapeDxQ + invJacXY .* shapeDyQ;
            shapeDyQnew = invJacXY .* shapeDxQ + invJacYY .* shapeDyQ;
            shapeDxQ = shapeDxQnew;
            shapeDyQ = shapeDyQnew;
            
            % Elemental matrices in PML
            weightsQ = weightsTriQ .* detJdxdu;
            matMpml = transpose(shapeOrQ) * (weightsQ .* shapeOrQ);
            matKpml = transpose(shapeDxQ) * (weightsQ .* shapeDxQ) + transpose(shapeDyQ) * (weightsQ .* shapeDyQ);
            matApml = (1/rho) * matKpml - (k^2/rho) * matMpml;
            
            % Elemental RHS vector in PML
            rhsProj = matMel\rhsPel;
            rhsPel = matApml * rhsProj - matAel * rhsProj;
            matAel = matApml;
        end
    end
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    rhsA(dof) = rhsA(dof) + rhsPel;
    
    iStart = (tri-1) * dofm.numDofPerTRI + 1;
    iEnd = tri * dofm.numDofPerTRI;
    matIv(iStart:iEnd,:) = dof' * ones(1,dofm.numDofPerTRI);
    matJv(iStart:iEnd,:) = ones(dofm.numDofPerTRI,1) * dof;
    matAv(iStart:iEnd,:) = matAel;
    % matMv(iStart:iEnd,:) = matMel;
    % matPv(iStart:iEnd,:) = matKel + k^2*matMel;
    
end
matA = sparse(matIv,matJv,matAv);
% matM = sparse(matIv,matJv,matMv);
% matP = sparse(matIv,matJv,matPv);
close(myWaitbar)
toc

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

tic
dofDIR0 = [];
dofDIR = [];
myWaitbar = waitbar(0,'   Surface terms...');
for edgBnd=1:mesh.numEdgBnd
    waitbar(edgBnd/mesh.numEdgBnd,myWaitbar,['   Surface terms... (' int2str(100*edgBnd/mesh.numEdgBnd) '%)']);
    edg = mesh.listEdgBnd(edgBnd);
    dof = dofm.locToGloBND(edgBnd,:);
    
    % Normal
    tri = mesh.mapEdgToTri(edg,1);
    fac = mesh.mapEdgToFac(edg,1);
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    normal = getNormalTRI(V1,V2,V3);
    normal = normal(fac,:);
    
    % Mapping
    ver = mesh.mapEdgToVer(edg,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    [xQ, yQ] = locToGloLIN(uLinQ, V1, V2);
    Jdxdu = norm(V2-V1) * 0.5;  % [ dx/du ]
    
    % Solution function
    [solQ, solDxQ, solDyQ] = mySourceSurface(xQ, yQ);
    dirQ = solQ;
    neuQ = normal(1)*solDxQ + normal(2)*solDyQ;
    
    % Orientation
    orientation = ones(dofm.numDofPerLIN,1);
    if(ver(1) > ver(2))
        orientation(3:dofm.numDofPerLIN) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    orientation = sparse(1:dofm.numDofPerLIN, 1:dofm.numDofPerLIN, orientation);
    
    % Shape function with orientation
    shapeOrQ = shapeLinQ * orientation;
    
    % Elemental matrices/vectors
    weightsQ = weightsLinQ * Jdxdu;
    matMel = shapeOrQ' * (weightsQ .* shapeOrQ);
    rhsDel = shapeOrQ' * (weightsQ .* dirQ);
    rhsNel = shapeOrQ' * (weightsQ .* neuQ);
    
    % Coefficients on the element
    k = kArray(tri);
    rho = rhoArray(tri);
    
    % Boundary condition
    switch edgTagToBC(mesh.tagEdgBnd(edgBnd))
        case 'DIR0'
            dofDIR0 = [dofDIR0 ; dof(:)];
        case 'DIR'
            dofDIR = [dofDIR ; dof(:)];
        case 'NEU0'
        case 'NEU'
            rhsA(dof) = rhsA(dof) + rhsNel;
        case 'ABC'
            matA(dof,dof) = matA(dof,dof) - 1i*k * matMel/rho;
        case 'ROB'
            matA(dof,dof) = matA(dof,dof) - 1i*k * matMel/rho;
            rhsA(dof) = rhsA(dof) + (rhsNel - 1i*k * rhsDel)/rho;
        otherwise
            error('BAD BOUNDARY CONDITION.');
    end
end
close(myWaitbar)

if(~isempty(dofDIR0))
    dofDIR0 = unique(dofDIR0);
    rhsA(dofDIR0) = 0;
    matA(dofDIR0,dofDIR0) = eye(size(dofDIR0,1),size(dofDIR0,1));
end
if(~isempty(dofDIR))
    dofDIR = unique(dofDIR);
    solP = computeSolProjL2_2D_CG(mesh, dofm);
    rhsA = rhsA - matA(:,dofDIR)*solP(dofDIR);
    rhsA(dofDIR) = solP(dofDIR);
    matA(dofDIR,dofDIR) = eye(size(dofDIR,1),size(dofDIR,1));
end
if(~isempty(pntSouTag))
    Isou = mesh.mapPntToVer(mesh.tagPntFile == pntSouTag);
    rhsA(Isou) = rhsA(Isou) + pntSouVal;
end
toc

% -------------------------------------------------------------------------
% Solve system
% -------------------------------------------------------------------------

% Matrix partition
numDofTRIred = mesh.numVer * dofm.numDofPerVer + mesh.numEdg * dofm.numDofPerEdg;
dofG = 1:numDofTRIred;                   % DOFs for reduced system
dofI = (numDofTRIred+1):dofm.numDofTRI;  % DOFs eliminated by static condensation
sysA.matII = matA(dofI,dofI);
sysA.matIG = matA(dofI,dofG);
sysA.matGI = matA(dofG,dofI);
sysA.matGG = matA(dofG,dofG);
sysA.rhsI = rhsA(dofI);
sysA.rhsG = rhsA(dofG);

% Full system
sysA.matA = matA;
sysA.rhsA = rhsA;

% Reduced system
% sysA.matS = sysA.matGG - sysA.matGI*(sysA.matII\sysA.matIG);
% sysA.rhsS = sysA.rhsG - sysA.matGI*(sysA.matII\sysA.rhsI);

% Preconditionning
if (PREC == 1)
    sysA.matP = 1;
    % sysA.matP = matM;
    % sysA.matP = matP;
else
    sysA.matP = 1;
end

% Compute solution
% solG = sysA.matS\sysA.rhsS;
% solI = sysA.matII\(sysA.rhsI-sysA.matIG*solG);
% solA = [ solG ; solI ];
tic
solA = sysA.matA\sysA.rhsA;
toc

end