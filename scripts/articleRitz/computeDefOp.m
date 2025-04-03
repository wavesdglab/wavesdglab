% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

function [Pdef,Qdef,Q] = computeDefOp(nbEigVec, eigenvec, A)

    if nbEigVec == 0
        Pdef = eye(size(A,1));
        Qdef = eye(size(A,1));
        Q = zeros(size(A,1));
        return;
    end

    Z = zeros(size(A,1),nbEigVec);
    Z(:,1:nbEigVec) = eigenvec;
    
    Z = orth(Z);
    
    Zt = Z';
    
    E = Zt*A*Z;
    
    temp = E\Zt;
    
    Q = Z*temp;
    
    Pdef = eye(size(A,1)) - A*Q;

    Qdef = eye(size(A,1)) - Q*A;

end