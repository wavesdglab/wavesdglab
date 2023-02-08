% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Assumption: x in [-1,1]

function val = functionsBernsteinDer(x,N)

x = x(:);
P = N-1;

val = zeros(size(x(:),1),N);
for i=0:P
    coef = factorial(P) / (factorial(i)*factorial(P-i));
    switch i
        case 0
            val(:,i+1) = coef * (P-i) * (-0.5) * (0.5-x/2).^(P-i-1);
        case P
            val(:,i+1) = coef * i * (0.5) * (0.5+x/2).^(i-1);
        otherwise
            val(:,i+1) = coef * (P-i) * (-0.5) * (0.5+x/2).^i .* (0.5-x/2).^(P-i-1) ...
                       + coef * i * (0.5) * (0.5+x/2).^(i-1) .* (0.5-x/2).^(P-i);
    end
end

valOut(:,1) = val(:,1);
valOut(:,2) = val(:,N);
valOut(:,3:N) = val(:,2:N-1);
val = valOut;

end