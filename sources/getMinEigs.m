% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

% This function returns the i smallest (i.e. close to 0) eigenvalues of A 

function [LambdaMin] = getMinEigs(A, i)

[~, eigenval] = eigs(A,size(A,1));
eigenval = diag(eigenval);

LambdaMin = zeros(1,i);

for j = 1:i
    % k = find(abs(eigenval) == min(abs(eigenval))) such that eigenval not in LambdaMin
    k = find(abs(eigenval) == min(abs(eigenval(~ismember(eigenval,LambdaMin)))));
    LambdaMin(j) = eigenval(k);
end

end