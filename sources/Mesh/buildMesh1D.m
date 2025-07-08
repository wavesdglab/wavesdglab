% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function mesh = buildMesh1D(Coord0, Coord1, numE)

mesh.numV   = numE+1;
mesh.numE   = numE;
mesh.coordV = linspace(Coord0, Coord1, mesh.numV);
mesh.listE  = [(1:mesh.numV-1)' (2:mesh.numV)'];

end