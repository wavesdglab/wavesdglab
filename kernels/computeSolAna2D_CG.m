function solAna = computeSolAna2D_CG(mesh)

x = mesh.coord(:,1);
y = mesh.coord(:,2);
solAna = mySol(x,y);

end