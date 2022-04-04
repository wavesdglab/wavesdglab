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

% Print interpolation scheme

fprintf(file,'$InterpolationScheme\n');
fprintf(file,'"INTERPOLATION_SCHEME"\n');
fprintf(file,'1\n');
fprintf(file,'3\n');        % 3 is for TRI
fprintf(file,'2\n');        % 2 is coded as-is in gmsh
fprintf(file,'3 3\n');      % size of coefficient matrix
fprintf(file,'1 -1 -1\n');
fprintf(file,'0 1 0\n');
fprintf(file,'0 0 1\n');
fprintf(file,'3 2\n');      % size of monomial matrix
fprintf(file,'0 0\n');
fprintf(file,'1 0\n');
fprintf(file,'0 1\n');

% 

[uQ, vQ, weights] = quadratureGaussTRI(2*dofm.degree);
valShape = functionsShapeTRI(uQ,vQ,dofm.degree);
uR = (uQ+1)/2;
vR = (vQ+1)/2;

% --- Print the coefficient matrix
numDofPerTRI = (dofm.degree+1)*(dofm.degree+2)/2;
fprintf(file,[int2str(numDofPerTRI) ' ' int2str(numDofPerTRI) '\n']);
%disp([int2str(numDofPerTRI) ' 2\n']);

vdm = zeros(numDofPerTRI,numDofPerTRI);
n1=0;
for i1=0:dofm.degree
    for j1=0:(dofm.degree-i1)
        n1=n1+1;
        n2=0;
        for i2=0:dofm.degree
            for j2=0:(dofm.degree-i2)
                n2=n2+1;
                vdm(n1,n2) = weights(:)' * ((uQ.^i1.*vQ.^j1) .* conj(uQ.^i2.*vQ.^j2));
            end
        end
    end
end

rhs = zeros(numDofPerTRI,1);
n1=0;
for i1=0:dofm.degree
    for j1=0:(dofm.degree-i1)
        n1=n1+1;
        for n2=1:numDofPerTRI
            rhs(n1,n2) = weights(:)' * ((uQ.^i1.*vQ.^j1) .* conj(valShape(:,n2)));
        end
    end
end

val = vdm\rhs;

for i=1:numDofPerTRI
    for j=1:numDofPerTRI
        fprintf(file,[num2str(val(i,j)) ' ']);
    end
    fprintf(file,'\n');
end


% --- Print the monomial matrix
fprintf(file,[int2str(numDofPerTRI) ' 2\n']);
%disp([int2str(numDofPerTRI) ' 2\n']);
for i=0:dofm.degree
    for j=0:(dofm.degree-i)
        fprintf(file,[int2str(i) ' ' int2str(j) '\n']);
        %disp([int2str(i) ' ' int2str(j) '\n']);
    end
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