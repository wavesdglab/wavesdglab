% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function [normErr, normSol] = computeNormError2D_DG(mesh, dofm, vecSol, vecRef)

global LdomX LdomY Rdom

% Quadrature
degreeQ = 3*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);

% Shape functions (f, dfdu, dfdv)
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);

normSolU2 = 0;
normErrU2 = 0;
normSolV2 = 0;
normErrV2 = 0;
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
    dofVx = 1*dofm.numDofTRI + dofm.locToGloTRI(tri,:);
    dofVy = 2*dofm.numDofTRI + dofm.locToGloTRI(tri,:);
    
    % Approximate solution (and derivatives)
    solQ = shapeOrQ * vecSol(dofU);
    solVxQ = shapeOrQ * vecSol(dofVx);
    solVyQ = shapeOrQ * vecSol(dofVy);
    
    % Reference solution (and derivatives)
    if (exist('vecRef','var'))
        refQ = shapeOrQ * vecRef(dofU);
        refVxQ = shapeOrQ * vecRef(dofVx);
        refVyQ = shapeOrQ * vecRef(dofVy);
    else
        refQ = mySol(xQ, yQ);
        [refQ, ~, ~, ~, refVxQ, refVyQ] = mySol(xQ, yQ);
    end
    
    % Error fields
    errQ = solQ(:) - refQ(:);
    errVxQ = solVxQ(:) - refVxQ(:);
    errVyQ = solVyQ(:) - refVyQ(:);
    
    % Error values
    normSolU2 = normSolU2 + weights' * (refQ .* conj(refQ)) * detJdxdu;
    normErrU2 = normErrU2 + weights' * (errQ .* conj(errQ)) * detJdxdu;
    normSolV2 = normSolV2 + weights' * (refVxQ .* conj(refVxQ) + refVyQ .* conj(refVyQ)) * detJdxdu;
    normErrV2 = normErrV2 + weights' * (errVxQ .* conj(errVxQ) + errVyQ .* conj(errVyQ)) * detJdxdu;
    
end

normSol = sqrt(normSolU2 + normSolV2);
normErr = sqrt(normErrU2 + normErrV2)/normSol;

end