% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

% This function compute the number of the nbEigVec closest eigenvectors to k for the Laplacian problem on a square domain with homogeneous Dirichlet boundary conditions

function indices = computeCloseEigVec_cavity(nb, k)
    
    limit = 100;
    [M, N] = meshgrid(1:limit, 1:limit);
    
    
    diff = abs(M.^2 + N.^2 - k^2/pi^2);
    
    
    [~, sorted_indices] = sort(diff(:));
    
    
    indices = sorted_indices(1:nb+1);

    
    [m, n] = ind2sub(size(diff), indices);

    
    if diff(m(nb), n(nb)) == diff(m(nb+1), n(nb+1))
        indices = [m(1:nb+1), n(1:nb+1)];
    else
        indices = [m(1:nb), n(1:nb)];
    end
end