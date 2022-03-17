function mesh = buildMesh1D(Coord0, Coord1, numV)

mesh.numV   = numV;
mesh.numE   = numV-1;
mesh.coordV = linspace(Coord0, Coord1, numV);
mesh.listE  = [(1:numV-1)' (2:numV)'];

end