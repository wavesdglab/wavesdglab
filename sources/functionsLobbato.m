% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Assumption: x in [-1,1]

function val = functionsLobbato(x,degree)

legendreInt = functionsLegendreInt(x,degree);
x = x(:);

N = degree+1;
val = zeros(size(x(:),1),N);

% nodal modes
val(:,1) = (1-x)/2;
val(:,2) = (1+x)/2;

% edge modes
for n=2:(N-1)
    val(:,n+1) = sqrt(n - 1/2) * legendreInt(:,n);
end

end