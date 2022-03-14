function solAna = computeSolAnaDG(mesh)

x = mesh.coord(:,1);
y = mesh.coord(:,2);
x = x(mesh.mapTriToVer)';
y = y(mesh.mapTriToVer)';
x = x(:);
y = y(:);
solAna = mySol(x,y);

end