function [normErr, normSol] = computeNormError2D_DG(mesh, dofm, vecSol, vecRef)

global LdomX LdomY Rdom Options

% Quadrature
degreeQ = 3*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);

% Shape functions (f, dfdu, dfdv)
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);

normSolU2 = 0;
normErrU2 = 0;
if strcmp(Options.Error,'Energy')
    normSolV2 = 0;
    normErrV2 = 0;
end
for tri=1:mesh.numTri

    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uQ, vQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
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

    dofU = 0*dofm.numDofTRI + dofm.locToGloTRI(tri,:);
    if strcmp(Options.Error,'Energy')
        dofVx = 1*dofm.numDofTRI + dofm.locToGloTRI(tri,:);
        dofVy = 2*dofm.numDofTRI + dofm.locToGloTRI(tri,:);
    end

    % Approximate solution (and derivatives)
    solQ = shapeOrQ * vecSol(dofU);
    if strcmp(Options.Error,'Energy')
        solVxQ = shapeOrQ * vecSol(dofVx);
        solVyQ = shapeOrQ * vecSol(dofVy);
    end

    % Reference solution (and derivatives)
    if (exist('vecRef','var') && ~isempty(vecRef))
        refQ = shapeOrQ * vecRef(dofU);
        if strcmp(Options.Error,'Energy')
            refVxQ = shapeOrQ * vecRef(dofVx);
            refVyQ = shapeOrQ * vecRef(dofVy);
        end
    else
        [refQ, ~, ~, refVxQ, refVyQ] = mySol(xQ, yQ);
    end

    % Error fields
    errQ = solQ(:) - refQ(:);
    if strcmp(Options.Error,'Energy')
        errVxQ = solVxQ(:) - refVxQ(:);
        errVyQ = solVyQ(:) - refVyQ(:);
    end

    % Error values
    normSolU2 = normSolU2 + weights' * (refQ .* conj(refQ)) * detJdxdu;
    normErrU2 = normErrU2 + weights' * (errQ .* conj(errQ)) * detJdxdu;
    if strcmp(Options.Error,'Energy')
        normSolV2 = normSolV2 + weights' * (refVxQ .* conj(refVxQ) + refVyQ .* conj(refVyQ)) * detJdxdu;
        normErrV2 = normErrV2 + weights' * (errVxQ .* conj(errVxQ) + errVyQ .* conj(errVyQ)) * detJdxdu;
    end

end

switch Options.Error
    case 'Energy'
        normSol = sqrt(normSolU2 + normSolV2);
        normErr = sqrt(normErrU2 + normErrV2)/normSol;
    otherwise
        normSol = sqrt(normSolU2);
        normErr = sqrt(normErrU2)/normSol;
end

end
