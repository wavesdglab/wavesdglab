% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

% This function compute the number of the nbEigVec closest eigenvectors to k for the Laplacian problem on a square domain with homogeneous Dirichlet boundary conditions

function indices = computeCloseEigVec_openCavity(nb, k)
    
    limit = 100;
    N = 1:limit;
    
    
    diff = abs((N./0.4).^2 - k^2/pi^2);
    
    
    [~, sorted_indices] = sort(diff(:));
    
    
    indices = sorted_indices(1:nb+1);

    
    n = ind2sub(size(diff), indices);

    if diff(n(nb)) == diff(n(nb+1))
        indices = n(1:nb+1);
    else
        indices = n(1:nb);
    end
end