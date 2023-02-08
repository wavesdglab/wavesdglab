% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Assumption: x in [-1,1]

function val = functionsShape1D(x,degree)

val = functionsLobbato(x,degree);
%val = functionsBernstein(x,degree);

end
