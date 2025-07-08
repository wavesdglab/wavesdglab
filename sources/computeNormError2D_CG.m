% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Authors: Axel Modave, Timothée Raynaud

function [error, norm] = computeNormError2D_CG(mesh, dofm, vecSol, vecRef)

global LdomX LdomY Rdom Options

% Quadrature
degreeQ = 3*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);

% Shape functions (f, dfdu, dfdv)
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);
[shapeDuQ, shapeDvQ] = functionsShapeDerTRI(uQ, vQ, dofm.degree);

normL2sol2 = 0;
errorL2sol2 = 0;
if strcmp(Options.Error,'H1')
    normL2der2 = 0;
    errorL2der2 = 0;
end
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

    if(~isempty(LdomX) && ~isempty(LdomY))
        ver = mesh.mapTriToVer(tri,:);
        VX = mesh.coord(ver,1);
        VY = mesh.coord(ver,2);
        if ((mean(abs(VX)) >= LdomX) || (mean(abs(VY)) >= LdomY))
            continue;
        end
    end
    if(~isempty(Rdom))
        ver = mesh.mapTriToVer(tri,:);
        VX = mesh.coord(ver,1);
        VY = mesh.coord(ver,2);
        VZ = VX + 1i*VY;
        VR = abs(VZ);
        if (mean(VR) >= Rdom)
            continue;
        end
    end

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
    if strcmp(Options.Error,'H1')
        shapeDxQ = (shapeDuQ * Jdudx(1,1) + shapeDvQ * Jdudx(2,1)) * orientation;
        shapeDyQ = (shapeDuQ * Jdudx(1,2) + shapeDvQ * Jdudx(2,2)) * orientation;
    end

    % Approximate solution (and derivatives)
    solQ   = shapeOrQ * vecSol(dofm.locToGloTRI(tri,:));
    if strcmp(Options.Error,'H1')
        solDxQ = shapeDxQ * vecSol(dofm.locToGloTRI(tri,:));
        solDyQ = shapeDyQ * vecSol(dofm.locToGloTRI(tri,:));
    end

    % Reference solution (and derivatives)
    if (exist('vecRef','var') && ~isempty(vecRef))
        refQ = shapeOrQ * vecRef(dofm.locToGloTRI(tri,:));
        if strcmp(Options.Error,'H1')
            refDxQ = shapeDxQ * vecRef(dofm.locToGloTRI(tri,:));
            refDyQ = shapeDyQ * vecRef(dofm.locToGloTRI(tri,:));
        end
    else
        refQ = mySol(xQ, yQ);
        if strcmp(Options.Error,'H1')
            [refQ, refDxQ, refDyQ] = mySol(xQ, yQ);
        end
    end

    % Error fields
    errQ = solQ(:) - refQ(:);
    if strcmp(Options.Error,'H1')
        errDxQ = solDxQ(:) - refDxQ(:);
        errDyQ = solDyQ(:) - refDyQ(:);
    end

    % Error values
    normL2sol2 = normL2sol2 + weights' * (refQ .* conj(refQ)) * detJdxdu;
    errorL2sol2 = errorL2sol2 + weights' * (errQ .* conj(errQ)) * detJdxdu;
    if strcmp(Options.Error,'H1')
        normL2der2 = normL2der2 + weights' * (refDxQ .* conj(refDxQ) + refDyQ .* conj(refDyQ)) * detJdxdu;
        errorL2der2 = errorL2der2 + weights' * (errDxQ .* conj(errDxQ) + errDyQ .* conj(errDyQ)) * detJdxdu;
    end
end

switch Options.Error
    case 'H1'
        norm = sqrt(normL2sol2 + normL2der2);
        error = sqrt(errorL2sol2 + errorL2der2)/norm;
    otherwise % L2
        norm = sqrt(normL2sol2);
        error = sqrt(errorL2sol2)/norm;
end
end