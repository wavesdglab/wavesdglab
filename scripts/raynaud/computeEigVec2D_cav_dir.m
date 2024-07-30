% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

% This function computes the nbEigVec closest eigenvectors to k for the Laplacian problem on a rectagular domain with homogeneous Dirichlet boundary conditions

function [eigenvec,nbEigVec] = computeEigVec2D_cav_dir(mesh, dofm, nbEigVec, varargin)

global k


if varargin{1} == "firsteigvec"
    mn = smallest_sum_of_squares(nbEigVec);
elseif varargin{1} == "closesteigvec"
    mn = closest_to_k(nbEigVec, k);
else
    error('Invalid input: varargin{1} must be "firsteigvec" or "closesteigvec"')
end

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
            n = mn(i);

            refQ(:, i) = (sin((xQ+0.75)*pi/1.3).*sin(n*(yQ-0.2)*pi/0.4));
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
%     filename = 'output/eigenvec' + string(i) + '.pos';
%     fieldname = 'eigenvec' + string(i);
%     writeField2D(dofm, mesh, eigenvec(:,i), filename, fieldname);
%     fprintf('field [%i] saved \n', i);
% end

end


function indices = closest_to_k(nb, k)
        % Generate a grid of integer indices up to a certain limit
        limit = 100;
        N = 1:limit;
        
        % Compute the difference between m^2 + n^2 and k^2/pi^2
        diff = abs((N./0.4).^2 - k^2/pi^2);
        
        % Flatten the matrices and sort the indices based on the difference
        [~, sorted_indices] = sort(diff(:));
        
        % Extract the nb+1 smallest indices
        indices = sorted_indices(1:nb+1);

        % Convert linear indices to subscripts
        n = ind2sub(size(diff), indices);

        if diff(n(nb)) == diff(n(nb+1))
            indices = n(1:nb+1);
        else
            indices = n(1:nb);
        end
end