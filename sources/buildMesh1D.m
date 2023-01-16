function mesh = buildMesh1D(Coord0, Coord1, numE)

mesh.numV   = numE+1;
mesh.numE   = numE;
mesh.coordV = linspace(Coord0, Coord1, mesh.numV);
mesh.listE  = [(1:mesh.numV-1)' (2:mesh.numV)'];

end