function writeField(dofm, mesh, field, nameFile, nameField)

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

% Compute/Print interpolation scheme

fprintf(file,'$InterpolationScheme\n');
fprintf(file,'"INTERPOLATION_SCHEME"\n'); % name(string)
fprintf(file,'1\n');                      % numElementTopologies(ASCII int)
fprintf(file,'3\n');                      % elementTopology - 3 is for TRI
fprintf(file,'2\n');                      % numInterpolationMatrices(ASCII int) - 2 is coded as-is in gmsh

% --- Compute matrices
[uQ, vQ, weights] = quadratureGaussTRI(4*dofm.degree);  % Quadrature
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);        % Shape functions
modalQ = zeros(size(uQ,1), dofm.numDofPerTRI);          % Modal functions
matExp = zeros(dofm.numDofPerTRI,2);
n = 1;
for n1 = 0:dofm.degree
    for n2 = 0:dofm.degree-n1
        matExp(n,1) = n1;
        matExp(n,2) = n2;
        modalQ(:,n) = (0.5*uQ+0.5).^n1 .* (0.5*vQ+0.5).^n2;
        n=n+1;
    end
end

for i=1:dofm.numDofPerTRI
    for j=1:dofm.numDofPerTRI
        matP(i,j) = weights' * (modalQ(:,i) .* modalQ(:,j));
        rhsP(i,j) = weights' * (modalQ(:,i) .* shapeQ(:,j));
    end
end
matCoef = matP\rhsP;
matCoef = matCoef';

% --- Print matrix of coefficients
fprintf(file,'%g %g\n', dofm.numDofPerTRI, dofm.numDofPerTRI);
for i=1:dofm.numDofPerTRI
    for j=1:dofm.numDofPerTRI
        fprintf(file,'%g ', matCoef(i,j));
    end
    fprintf(file,'\n');
end

% --- Print matrix of exponants
fprintf(file,'%g %g\n', dofm.numDofPerTRI, 2);
for i=1:dofm.numDofPerTRI
    for j=1:2
        fprintf(file,'%g ', matExp(i,j));
    end
    fprintf(file,'\n');
end

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
for tri=1:mesh.numTri
    fprintf(file,'%i %i ', tri, dofm.numDofPerTRI);
    for n=1:dofm.numDofPerTRI
        fprintf(file,'%f ', real(field(dofm.locToGloTRI(tri,n))));
    end
    fprintf(file,'\n');
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
fprintf(file,'1\n'); % IMAG PART
fprintf(file,'1\n');
fprintf(file,'%i\n',mesh.numTri);
fprintf(file,'0\n');
for tri=1:mesh.numTri
    fprintf(file,'%i %i ', tri, dofm.numDofPerTRI);
    for n=1:dofm.numDofPerTRI
        fprintf(file,'%f ', imag(field(dofm.locToGloTRI(tri,n))));
    end
    fprintf(file,'\n');
end
fprintf(file,'$EndElementNodeData\n');

% Close file

fclose(file);

end