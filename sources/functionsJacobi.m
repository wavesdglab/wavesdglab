% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Assumption: x in [-1,1]

function val = functionsJacobi(x,alpha,beta,N)

x = x(:);

val = zeros(size(x,1),N);
val(:,1) = 1;                                  % order 0
val(:,2) = 0.5*(alpha-beta+(alpha+beta+2)*x);  % order 1
for n=3:N
    a1 = 2*(n-1)*(n+alpha+beta-1)*(2*n+alpha+beta-4);
    a2 = (2*n+alpha+beta-3)*(alpha^2-beta^2);
    a3 = (2*n+alpha+beta-3)*(2*n+alpha+beta-4)*(2*n+alpha+beta-2);
    a4 = 2*(n+alpha-2)*(n+beta-2)*(2*n+alpha+beta-2);
    val(:,n) = ((a2+a3*x).*val(:,n-1) - a4*val(:,n-2)) / a1;  % order (n-1)
end

end