% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

function [P,Q] = computeDefOp(nbEigVec, eigenvec, A)


    Z = zeros(size(A,1),nbEigVec);
    Z(:,1:nbEigVec) = eigenvec;
    
    Z = orth(Z);
    
    Zt = Z';
    
    E = Zt*A*Z;
    
    temp = E\Zt;
    
    Q = Z*temp;
    
    P = eye(size(A,1)) - A*Q;

end