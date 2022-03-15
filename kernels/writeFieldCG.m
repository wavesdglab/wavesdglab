function writeFieldCG(mesh, field, nameFile, nameField)

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

% Print interpolation scheme

fprintf(file,'$InterpolationScheme\n');
fprintf(file,'"INTERPOLATION_SCHEME"\n');
fprintf(file,'1\n');
fprintf(file,'3\n');
fprintf(file,'2\n');
fprintf(file,'3 3\n');
fprintf(file,'1 -1 -1\n');
fprintf(file,'0 1 0\n');
fprintf(file,'0 0 1\n');
fprintf(file,'3 2\n');
fprintf(file,'0 0\n');
fprintf(file,'1 0\n');
fprintf(file,'0 1\n');
fprintf(file,'$EndInterpolationScheme\n');

% Print element node data (real part)

fprintf(file,'$ElementNodeData\n');
fprintf(file,'2\n');
fprintf(file,'"%s"\n',nameField);
fprintf(file,'"INTERPOLATION_SCHEME"\n');
fprintf(file,'1\n');
fprintf(file,'0\n');
fprintf(file,'4\n');
fprintf(file,'0\n'); % REAL PART
fprintf(file,'1\n');
fprintf(file,'%i\n',mesh.numTri);
fprintf(file,'0\n');
for i=1:mesh.numTri
    fprintf(file,'%i 3 %f %f %f\n', i, ...
        real(field(mesh.mapTriToVer(i,1))), ...
        real(field(mesh.mapTriToVer(i,2))), ...
        real(field(mesh.mapTriToVer(i,3))));
end
fprintf(file,'$EndElementNodeData\n');

% Print element node data (imaginary part)

fprintf(file,'$ElementNodeData\n');
fprintf(file,'2\n');
fprintf(file,'"%s"\n',nameField);
fprintf(file,'"INTERPOLATION_SCHEME"\n');
fprintf(file,'1\n');
fprintf(file,'0\n');
fprintf(file,'4\n');
fprintf(file,'1\n'); % REAL PART
fprintf(file,'1\n');
fprintf(file,'%i\n',mesh.numTri);
fprintf(file,'0\n');
for i=1:mesh.numTri
    fprintf(file,'%i 3 %f %f %f\n', i, ...
        imag(field(mesh.mapTriToVer(i,1))), ...
        imag(field(mesh.mapTriToVer(i,2))), ...
        imag(field(mesh.mapTriToVer(i,3))));
end
fprintf(file,'$EndElementNodeData\n');

% Close file

fclose(file);

end