function solAna = computeSolAnaCG(mesh)

x = mesh.coord(:,1);
y = mesh.coord(:,2);
solAna = mySol(x,y);

end