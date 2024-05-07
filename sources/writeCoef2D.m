% Copyright (C) 2024, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function writeCoef2D(mesh, field, nameFile, nameField)

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

% Print element node data (real part)

fprintf(file,'$ElementData\n');
fprintf(file,'2\n');
fprintf(file,'"%s"\n',nameField);
fprintf(file,'""\n');
fprintf(file,'1\n');
fprintf(file,'0\n');
fprintf(file,'4\n');
fprintf(file,'0\n'); % timeStep
fprintf(file,'1\n'); % numComp
fprintf(file,'%i\n',mesh.numTri); % numEnt
fprintf(file,'0\n'); % partition
for tri=1:mesh.numTri
    fprintf(file,'%i %f\n', tri, real(field(tri)));
end
fprintf(file,'$EndElementData\n');

% Close file

fclose(file);

end