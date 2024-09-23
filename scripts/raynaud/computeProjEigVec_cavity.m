% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

% This function computes the nbEigVec first eigenvectors or the nbEigVec closest eigenvectors to k for the Laplacian problem on a square domain with homogeneous Dirichlet boundary conditions

function [eigenvec,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec, varargin)

global k


if varargin{1} == "firsteigvec"
    mn = smallest_sum_of_squares(nbEigVec);
elseif varargin{1} == "closesteigvec"
    mn = closest_to_k(nbEigVec, k);
else
    error('Invalid input: varargin{1} must be "firsteigvec" or "closesteigvec"')
end

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
    for i=1:size(mn, 1)
        m = mn(i, 1);
        n = mn(i, 2);

        refQ(:, i) = (sin(m*xQ*pi).*sin(n*yQ*pi)) .* (16*m*n*pi^2 /(pi^2*(m^2 + n^2)-k^2));
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


    nbEigVec = size(mn, 1);
    
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

function indices = smallest_sum_of_squares(nb)
        % Generate a grid of integer indices up to a certain limit
        limit = 100;
        [M, N] = meshgrid(1:limit, 1:limit);
        
        % Compute the sum of squares m^2 + n^2
        sum_squares = M.^2 + N.^2;
        
        % Flatten the matrices and sort the indices based on sum of squares
        [~, sorted_indices] = sort(sum_squares(:));
        
        % Extract the nb+1 smallest indices
        nb_smallest_indices = sorted_indices(1:nb+1);
        
        % Convert linear indices to subscripts
        [m, n] = ind2sub(size(sum_squares), nb_smallest_indices);

        %% test if the nbth and the nb+1th smallest are equal
        if sum_squares(m(nb), n(nb)) == sum_squares(m(nb+1), n(nb+1))
            indices = [m(1:nb+1), n(1:nb+1)];
        else
            indices = [m(1:nb), n(1:nb)];
        end
end

function indices = closest_to_k(nb, k)
        % Generate a grid of integer indices up to a certain limit
        limit = 100;
        [M, N] = meshgrid(1:limit, 1:limit);
        
        % Compute the difference between m^2 + n^2 and k^2/pi^2
        diff = abs(M.^2 + N.^2 - k^2/pi^2);
        
        % Flatten the matrices and sort the indices based on the difference
        [~, sorted_indices] = sort(diff(:));
        
        % Extract the nb+1 smallest indices
        indices = sorted_indices(1:nb+1);

        % Convert linear indices to subscripts
        [m, n] = ind2sub(size(diff), indices);

        %% test if the nbth and the nb+1th smallest are equal
        if diff(m(nb), n(nb)) == diff(m(nb+1), n(nb+1))
            indices = [m(1:nb+1), n(1:nb+1)];
        else
            indices = [m(1:nb), n(1:nb)];
        end
end