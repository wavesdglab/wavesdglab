% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

% This function computes the nbEigVec closest eigenvectors to k for the Laplacian problem on a rectagular domain with homogeneous Neumann boundary conditions

function [eigenvec,nbEigVec] = computeProjEigVec_openCavity_NEU(mesh, dofm, nbEigVec, k)

mn = computeCloseEigVec_openCavity_NEU(nbEigVec, k);

nbEigVec = size(mn, 1);

% Quadrature and shape functions
degreeQ = 4*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);

% Build matrix and RHS vector
matP = sparse(dofm.numDofTRI, dofm.numDofTRI);
rhsP = zeros(dofm.numDofTRI, size(mn, 1));

for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uQ, vQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    detJdxdu = abs(det(Jdxdu));
    
    % Reference solution
    refQ = zeros(size(xQ,1), size(mn, 1));
    
    if (V1(1) >= -0.75 && V2(1) >= -0.75 && V3(1) >= -0.75 && V1(1) <= 0.55 && V2(1) <= 0.55 && V3(1) <= 0.55 && V1(2) >= -0.2 && V2(2) >= -0.2 && V3(2) >= -0.2 && V1(2) <= 0.2 && V2(2) <= 0.2 && V3(2) <= 0.2)
        for i=1:size(mn, 1)
            m = mn(i, 1);
            n = mn(i, 2);
            refQ(:, i) = (sin((m+1/2)*(xQ+0.75)*pi/1.3).*cos(n*(yQ-0.2)*pi/0.4));
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
    
    % Shape functions with orientation
    shapeOrQ = shapeQ * orientation;
    
    % Local matrix and RHS vector
    matPel = shapeOrQ' * weights * shapeOrQ * detJdxdu;
    rhsPel = shapeOrQ' * weights * refQ * detJdxdu;
    
    % Assembling
    dof = dofm.locToGloTRI(tri,:);
    matP(dof,dof) = matP(dof,dof) + matPel;
    rhsP(dof,:) = rhsP(dof,:) + rhsPel;
    
end

% Solution
eigenvec = matP\rhsP;


% for i=1:size(mn, 1)
%     m = mn(i, 1);
%     n = mn(i, 2);
%     filename = 'output/eigenvec' + string(m)+'*'+string(n) + '.pos';
%     fieldname = 'eigenvec' + string(m)+'*'+string(n);
%     writeField2D(dofm, mesh, eigenvec(:,i), filename, fieldname);
%     fprintf('field [%i] saved \n', i);
    
% end

end


% This function compute the number of the nbEigVec closest eigenvectors to
% k for the Laplacian problem on a rectangular domain with homogeneous
% Neumamnn BC

function indices = computeCloseEigVec_openCavity_NEU(nb, k)

limit = 100;
M = 0:limit-1;
N = 1:limit;
[M, N] = meshgrid(M, N);

quasi_resonances = (M+1/2).^2/((1.3)^2) + (N.^2)/(0.4^2);

diff = abs(quasi_resonances - k^2/pi^2);

[~, sorted_indices] = sort(diff(:));


indices = sorted_indices(1:nb);


[n, m] = ind2sub(size(diff), indices);

m = m - 1;

indices = [m(1:nb), n(1:nb)];
end
