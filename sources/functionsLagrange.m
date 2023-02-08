% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Assumption: x in [-1,1]

function val = functionsLagrange(x,nodes)

x = x(:);
nodes = nodes(:);
N = size(nodes,1);

val = ones(size(x(:),1),N);
for n=1:N
    for i=1:N
        if i~=n
            val(:,n) = val(:,n) .* (x - nodes(i)) / (nodes(n) - nodes(i));
        end
    end
end

end