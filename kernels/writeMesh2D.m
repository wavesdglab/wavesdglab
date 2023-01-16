function writeMesh2D(mesh, nameFile)

% Open file

file = fopen(nameFile,'w');

% Print mesh format

fprintf(file,'$MeshFormat\n');
fprintf(file,'4.1 0 8\n');
fprintf(file,'$EndMeshFormat\n');

% Print entities

fprintf(file,'$Entities\n');
fprintf(file,'0 0 1 0\n');
fprintf(file,'1 -0.99 0.99 -0.99 0.99 -0.99 0.99 1 888 0\n'); % 999 is surfaceTag ; 888 is physicalTag
fprintf(file,'$EndEntities\n');

% Print nodes

fprintf(file,'$Nodes\n');
fprintf(file,'1 %i 1 %i\n', mesh.numVer, mesh.numVer);
fprintf(file,'2 1 0 %i\n', mesh.numVer);
for i=1:mesh.numVer
    fprintf(file,'%i\n', i);
end
for i=1:mesh.numVer
    fprintf(file,'%f %f %f\n', mesh.coord(i,1), mesh.coord(i,2), 0.);
end
fprintf(file,'$EndNodes\n');

% Print elements

fprintf(file,'$Elements\n');
fprintf(file,'1 %i 1 %i\n', mesh.numTri, mesh.numTri);
fprintf(file,'2 1 2 %i\n', mesh.numTri);
for i=1:mesh.numTri
    fprintf(file,'%i %i %i %i\n', i, mesh.mapTriToVer(i,1), mesh.mapTriToVer(i,2), mesh.mapTriToVer(i,3));
end
fprintf(file,'$EndElements\n');

% Close file

fclose(file);

end