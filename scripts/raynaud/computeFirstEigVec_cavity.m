% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

% This function compute the number of the nbEigVec first eigenvectors for the Laplacian problem on a square domain with homogeneous Dirichlet boundary conditions


function indices = computeFirstEigVec_cavity(nb)
    
    limit = 100;
    [M, N] = meshgrid(1:limit, 1:limit);
    
    
    sum_squares = M.^2 + N.^2;
    
    
    [~, sorted_indices] = sort(sum_squares(:));
    
    
    nb_smallest_indices = sorted_indices(1:nb+1);
    
    
    [m, n] = ind2sub(size(sum_squares), nb_smallest_indices);

    
    if sum_squares(m(nb), n(nb)) == sum_squares(m(nb+1), n(nb+1))
        indices = [m(1:nb+1), n(1:nb+1)];
    else
        indices = [m(1:nb), n(1:nb)];
    end
end
