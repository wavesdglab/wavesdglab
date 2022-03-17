function solAna = computeSolAna2D_DG(mesh)

x = mesh.coord(:,1);
y = mesh.coord(:,2);
x = x(mesh.mapTriToVer)';
y = y(mesh.mapTriToVer)';
x = x(:);
y = y(:);
solAna = mySol(x,y);

end