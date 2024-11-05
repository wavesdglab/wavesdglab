% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Authors: Axel Modave, Timothée Raynaud

function [solA, sysA] = computeSolNum2D_CG(mesh, dofm, PREC)

global k edgTagToBC
global LdomX LdomY LpmlX LpmlY Rdom Rpml

matA = sparse(dofm.numDofTRI,dofm.numDofTRI);
matM = sparse(dofm.numDofTRI,dofm.numDofTRI);
matShiftedLaplacian = sparse(dofm.numDofTRI,dofm.numDofTRI);
rhsA = zeros(dofm.numDofTRI, 1);

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
    
    % RHS function
    rhsQ = mySourceVolume(xQ, yQ);
    
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
        end
    end
    
    % Elemental matrices/vectors
    weightsQ = weightsTriQ .* detJdxdu;
    matMel = transpose(shapeOrQ) * (weightsQ .* shapeOrQ);
    matKel = transpose(shapeDxQ) * (weightsQ .* shapeDxQ) + transpose(shapeDyQ) * (weightsQ .* shapeDyQ);
    rhsPel = transpose(shapeOrQ) * (weightsQ .* rhsQ);
    
    % PML stretching (rectangular PML — alternative)
    % if(~isempty(LdomX) && ~isempty(LdomY))
    %     if ((mean(abs(xQ)) >= LdomX) || (mean(abs(yQ)) >= LdomY))
    %         sigmaPmlX = (LdomX <= abs(xQ))./(LdomX+LpmlX-abs(xQ));
    %         sigmaPmlY = (LdomY <= abs(yQ))./(LdomY+LpmlY-abs(yQ));
    %         gammaPmlX = ones(size(xQ)) - sigmaPmlX/(1i*k);
    %         gammaPmlY = ones(size(yQ)) - sigmaPmlY/(1i*k);
    %         coefPml = gammaPmlX.*gammaPmlY;
    %         tensPmlXX = gammaPmlY./gammaPmlX;
    %         tensPmlYY = gammaPmlX./gammaPmlY;
    %         matMel = transpose(shapeOrQ) * (weightsQ .* coefPml .* shapeOrQ);
    %         matKel = transpose(shapeDxQ) * (weightsQ .* tensPmlXX .* shapeDxQ) + transpose(shapeDyQ) * (weightsQ .* tensPmlYY .* shapeDyQ);
    %         rhsPel = transpose(shapeOrQ) * (weightsQ .* rhsQ);
    %     end
    % end
    
    % % PML stretching (circular PML — alternative)
    % if(~isempty(Rdom))
    %     rQ = sqrt(xQ.*xQ + yQ.*yQ);
    %     if (mean(rQ) >= Rdom)
    %         cosT = xQ./rQ;
    %         sinT = yQ./rQ;
    %         sigmaPml = 1./(Rpml-(rQ-Rdom));
    %         sigmaPmlInt = -log(1-(rQ-Rdom)/Rpml);
    %         gammaPmlR = ones(size(rQ)) - sigmaPml/(1i*k);
    %         gammaPmlT = ones(size(rQ)) - sigmaPmlInt/(1i*k)./rQ;
    %         coefPml = gammaPmlR.*gammaPmlT;
    %         tensPmlXX = (gammaPmlT./gammaPmlR) .* cosT.*cosT + (gammaPmlR./gammaPmlT) .* (sinT.*sinT);
    %         tensPmlXY = (gammaPmlT./gammaPmlR) .* cosT.*sinT - (gammaPmlR./gammaPmlT) .* (cosT.*sinT);
    %         tensPmlYY = (gammaPmlT./gammaPmlR) .* sinT.*sinT + (gammaPmlR./gammaPmlT) .* (cosT.*cosT);
    %         matMel = transpose(shapeOrQ) * (weightsQ .* coefPml .* shapeOrQ);
    %         matKel = transpose(shapeDxQ) * (weightsQ .* tensPmlXX .* shapeDxQ) ...
    %                + transpose(shapeDxQ) * (weightsQ .* tensPmlXY .* shapeDyQ) ...
    %                + transpose(shapeDyQ) * (weightsQ .* tensPmlXY .* shapeDxQ) ...
    %                + transpose(shapeDyQ) * (weightsQ .* tensPmlYY .* shapeDyQ);
    %         rhsPel = transpose(shapeOrQ) * (weightsQ .* rhsQ);
    %         
    %     end
    % end
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    matA(dof,dof) = matA(dof,dof) + matKel - k^2*matMel;
    matM(dof,dof) = matM(dof,dof) + matMel;
    rhsA(dof) = rhsA(dof) + rhsPel;
    matShiftedLaplacian(dof,dof) = matShiftedLaplacian(dof,dof) + matKel - (k*k + 1i*k)*matMel;
    
end

% -------------------------------------------------------------------------
% Surface terms
% -------------------------------------------------------------------------

dofDIR = [];
cacheDIR = zeros(dofm.numDofTRI);
for edgBnd=1:mesh.numEdgBnd
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
    
    % Boundary condition
    switch edgTagToBC(mesh.tagEdgBnd(edgBnd))
        case 'DIR0'
            dofDIR = [dofDIR ; dof];
            cacheDIR(dof) = zeros(size(dof,1),size(dof,2));
        case 'DIR'
            dofDIR = [dofDIR ; dof];
            cacheDIR(dof) = ones(size(dof,1),size(dof,2));
        case 'NEU0'
        case 'NEU'
            rhsA(dof) = rhsA(dof) + rhsNel;
        case 'ABC'
            matA(dof,dof) = matA(dof,dof) - 1i*k * matMel;
        case 'ROB'
            matA(dof,dof) = matA(dof,dof) - 1i*k * matMel;
            rhsA(dof) = rhsA(dof) + rhsNel - 1i*k * rhsDel;
        otherwise
            error('BAD BOUNDARY CONDITION.');
    end
end

if(~isempty(dofDIR))
    solP = computeSolProjL2_2D_CG(mesh, dofm) .* cacheDIR;
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
sysA.matS = sysA.matGG - sysA.matGI*(sysA.matII\sysA.matIG);
sysA.rhsS = sysA.rhsG - sysA.matGI*(sysA.matII\sysA.rhsI);

% Preconditionning
if (PREC == 1)
    % sysA.matP = matM;
    sysA.matP = matShiftedLaplacian;
else
    sysA.matP = 1;
end

% Compute solution
solG = sysA.matS\sysA.rhsS;
solI = sysA.matII\(sysA.rhsI-sysA.matIG*solG);
solA = [ solG ; solI ];
% solA = sysA.matA\sysA.rhsA;

end