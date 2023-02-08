% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Assumption: x in [-1,1]

function val = functionsBernstein(x,N)

x = x(:);
P = N-1;

val = zeros(size(x(:),1),N);
for i=0:P
    coef = factorial(P) / (factorial(i)*factorial(P-i));
    val(:,i+1) = coef * (0.5+x/2).^i .* (0.5-x/2).^(P-i); % order i
end

valOut(:,1) = val(:,1);
valOut(:,2) = val(:,N);
valOut(:,3:N) = val(:,2:N-1);
val = valOut;

end