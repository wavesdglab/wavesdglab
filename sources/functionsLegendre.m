% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Assumption: x in [-1,1]

function val = functionsLegendre(x,degree)

x = x(:);

N = degree+1;
val = zeros(size(x(:),1),N);
val(:,1) = 1;  % order 0
val(:,2) = x;  % order 1
for n=2:(N-1)
    val(:,n+1) = ((2*n-1)*x.*val(:,n) - (n-1)*val(:,n-1)) / n;  % order n
end
for n=1:N
    val(:,n) = val(:,n) * sqrt(n-0.5);  % order n
end

end