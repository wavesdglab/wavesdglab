% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Assumption: x in [-1,1]

function val = functionsShapeDer1D(x,degree)

val = functionsLobbatoDer(x,degree);
%val = functionsBernsteinDer(x,degree);

end