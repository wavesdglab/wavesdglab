% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

% This function computes the nbEigVec closest eigenvectors to k for the Laplacian problem on a rectagular domain with homogeneous Dirichlet boundary conditions

function [eigenvec,nbEigVec,rhsP,matP] = computeProjEigVec_openCavity_DIR(mesh, dofm, nbEigVec, k)

global edgTagToBC Options

if nbEigVec == 0
    eigenvec = [];
    return;
end

mn = computeCloseEigVec_openCavity_DIR(nbEigVec, k);

nbEigVec = size(mn, 1);

% Quadrature and shape functions
degreeQ = 2*dofm.degree;
[uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);
shapeQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);

% Build matrix and RHS vector
matP = sparse(dofm.numDofTRI, dofm.numDofTRI);
rhsP = zeros(dofm.numDofTRI, size(mn, 1));
for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uTriQ, vTriQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    detJdxdu = abs(det(Jdxdu));
    
    % Reference solution
    refQ = zeros(size(xQ,1), size(mn, 1));
    
    eps = 0.05;
    ymin = -0.2 - eps;
    ymax = 0.2 + eps;
    xmin = -0.75 - eps;
    xmax = 0.55 + eps;

    if (V1(1) >= xmin && V2(1) >= xmin && V3(1) >= xmin && V1(1) <= xmax && V2(1) <= xmax && V3(1) <= xmax && V1(2) >= ymin && V2(2) >= ymin && V3(2) >= ymin && V1(2) <= ymax && V2(2) <= ymax && V3(2) <= ymax)
        for i=1:size(mn, 1)
            m = mn(i, 1);
            n = mn(i, 2);
            refQ(:, i) = (sin(m*(xQ+0.75)*pi/1.3).*sin(n*(yQ-0.2)*pi/0.4));
        end
    end
    
    % Orientation
    orientation = ones(dofm.numDofPerTRI,1);
    if ~strcmp(Options.Basis,'Lagrange')
        if(ver(1) > ver(2))
            orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        if(ver(2) > ver(3))
            orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        if(ver(3) > ver(1))
            orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
    end
    orientation = sparse(1:dofm.numDofPerTRI, 1:dofm.numDofPerTRI, orientation);

    % Shape functions with orientation
    shapeOrQ = shapeQ * orientation;
    
    % Local matrix and RHS vector
    weightsQ = weightsTriQ .* detJdxdu;
    matPel = transpose(shapeOrQ) * (weightsQ .* shapeOrQ);
    rhsPel = transpose(shapeOrQ) * (weightsQ .* refQ);
    
    % Assembling
    dof = dofm.locToGloTRI(tri,:);
    matP(dof,dof) = matP(dof,dof) + matPel;
    rhsP(dof,:) = rhsP(dof,:) + rhsPel;

end


dofDIR = [];
cacheDIR = zeros(dofm.numDofTRI);
for edgBnd=1:mesh.numEdgBnd
    dof = dofm.locToGloBND(edgBnd,:);
    % Boundary condition
    switch edgTagToBC(mesh.tagEdgBnd(edgBnd))
        case 'DIR0'
            dofDIR = [dofDIR ; dof];
            cacheDIR(dof) = zeros(size(dof,1),size(dof,2));
    end
end
if(~isempty(dofDIR))
    dofDIR = unique(dofDIR);
    matP(:,dofDIR) = 0;
    matP(dofDIR,:) = 0;
    matP(dofDIR,dofDIR) = eye(size(dofDIR,1),size(dofDIR,1));
    rhsP(dofDIR,:) = 0;
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
% Dirichlet BC

function indices = computeCloseEigVec_openCavity_DIR(nb, k)

limit = 100;
M = 1:limit;
N = 1:limit;
[M, N] = meshgrid(M, N);

quasi_resonances = (M.^2)/1.3^2 + (N.^2)/0.4^2;

diff = abs(quasi_resonances - k^2/pi^2);

[~, sorted_indices] = sort(diff(:));


indices = sorted_indices(1:nb);


[n, m] = ind2sub(size(diff), indices);


indices = [m(1:nb), n(1:nb)];
end
