function writeField2D(dofm, mesh, field, nameFile, nameField)

global LdomX LdomY Rdom PML_HIDE

% Open file

file = fopen(nameFile,'w');

% Print mesh format

fprintf(file,'$MeshFormat\n');
fprintf(file,'4.1 0 8\n');
fprintf(file,'$EndMeshFormat\n');

% % Print entities
%
% fprintf(file,'$Entities\n');
% fprintf(file,'0 0 1 0\n');
% fprintf(file,'1 -0.99 0.99 -0.99 0.99 -0.99 0.99 1 888 0\n'); % 999 is surfaceTag ; 888 is physicalTag
% fprintf(file,'$EndEntities\n');
%
% % Print nodes
%
% fprintf(file,'$Nodes\n');
% fprintf(file,'1 %i 1 %i\n', mesh.numVer, mesh.numVer);
% fprintf(file,'2 1 0 %i\n', mesh.numVer);
% for i=1:mesh.numVer
%     fprintf(file,'%i\n', i);
% end
% for i=1:mesh.numVer
%     fprintf(file,'%f %f %f\n', mesh.coord(i,1), mesh.coord(i,2), 0.);
% end
% fprintf(file,'$EndNodes\n');
%
% % Print elements
%
% fprintf(file,'$Elements\n');
% fprintf(file,'1 %i 1 %i\n', mesh.numTri, mesh.numTri);
% fprintf(file,'2 1 2 %i\n', mesh.numTri);
% for tri=1:mesh.numTri
%     fprintf(file,'%i %i %i %i\n', mesh.tagTriFile(tri), mesh.mapTriToVer(tri,1), mesh.mapTriToVer(tri,2), mesh.mapTriToVer(tri,3));
% end
% fprintf(file,'$EndElements\n');

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

% If PML

numTriPml = 0;
if(PML_HIDE == 1)
    for tri=1:mesh.numTri
        ver = mesh.mapTriToVer(tri,:);
        if(~isempty(LdomX) && ~isempty(LdomY))
            VX = mesh.coord(ver,1);
            VY = mesh.coord(ver,2);
            if ((mean(abs(VX)) >= LdomX) || (mean(abs(VY)) >= LdomY))
                numTriPml = numTriPml+1;
            end
        end
        if(~isempty(Rdom))
            VX = mesh.coord(ver,1);
            VY = mesh.coord(ver,2);
            VZ = VX + 1i*VY;
            VR = abs(VZ);
            if (mean(VR) >= Rdom)
                numTriPml = numTriPml+1;
            end
        end
    end
end

% Print element node data [absolute value or real/imaginary values]

global WRITE_FIELD_ABSOLUTE
if(WRITE_FIELD_ABSOLUTE == 1)

    % Print element node data (absolute value)

    degreeQ = 2*dofm.degree;
    [uTriQ, vTriQ, weightsTriQ] = quadratureGaussTRI(degreeQ);
    shapeQ = functionsShapeTRI(uTriQ, vTriQ, dofm.degree);

    fprintf(file,'$ElementNodeData\n');
    fprintf(file,'2\n');
    fprintf(file,'"%s"\n',nameField);
    fprintf(file,'"INTERPOLATION_SCHEME"\n');
    fprintf(file,'1\n');
    fprintf(file,'0\n');
    fprintf(file,'4\n');
    fprintf(file,'2\n'); % IMAG PART
    fprintf(file,'1\n');
    fprintf(file,'%i\n',mesh.numTri-numTriPml);
    fprintf(file,'0\n');
    for tri=1:mesh.numTri
        ver = mesh.mapTriToVer(tri,:);

        if(PML_HIDE == 1)
            if(~isempty(LdomX) && ~isempty(LdomY))
                VX = mesh.coord(ver,1);
                VY = mesh.coord(ver,2);
                if ((mean(abs(VX)) >= LdomX) || (mean(abs(VY)) >= LdomY))
                    continue;
                end
            end
            if(~isempty(Rdom))
                VX = mesh.coord(ver,1);
                VY = mesh.coord(ver,2);
                VZ = VX + 1i*VY;
                VR = abs(VZ);
                if (mean(VR) >= Rdom)
                    continue;
                end
            end
        end

        % Orientation
        orientation = ones(dofm.numDofPerTRI,1);
        if(ver(1) > ver(2))
            orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        if(ver(2) > ver(3))
            orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        if(ver(3) > ver(1))
            orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        orientation = sparse(1:dofm.numDofPerTRI, 1:dofm.numDofPerTRI, orientation);

        % Absolute value of field at quadrature nodes
        shapeOrQ = shapeQ * orientation;  % Shape function
        solQ = shapeOrQ * field(dofm.locToGloTRI(tri,:));
        solQ = abs(solQ);

        % Absolute value of field
        V1 = mesh.coord(ver(1),:);
        V2 = mesh.coord(ver(2),:);
        V3 = mesh.coord(ver(3),:);
        Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
        detJdxdu = abs(det(Jdxdu));
        weightsQ = weightsTriQ .* detJdxdu;
        matPel = transpose(shapeOrQ) * (weightsQ .* shapeOrQ);
        rhsPel = transpose(shapeOrQ) * (weightsQ .* solQ);
        fieldTri = matPel\rhsPel;
        fieldTri = orientation*fieldTri;

        % Print
        fprintf(file,'%i %i ', mesh.tagTriFile(tri), dofm.numDofPerTRI);
        for n=1:dofm.numDofPerTRI
            fprintf(file,'%f ', fieldTri(n));
        end
        fprintf(file,'\n');
    end
    fprintf(file,'$EndElementNodeData\n');

else

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
    fprintf(file,'%i\n',mesh.numTri-numTriPml);
    fprintf(file,'0\n');
    for tri=1:mesh.numTri
        ver = mesh.mapTriToVer(tri,:);

        if(PML_HIDE == 1)
            if(~isempty(LdomX) && ~isempty(LdomY))
                VX = mesh.coord(ver,1);
                VY = mesh.coord(ver,2);
                if ((mean(abs(VX)) >= LdomX) || (mean(abs(VY)) >= LdomY))
                    continue;
                end
            end
            if(~isempty(Rdom))
                VX = mesh.coord(ver,1);
                VY = mesh.coord(ver,2);
                VZ = VX + 1i*VY;
                VR = abs(VZ);
                if (mean(VR) >= Rdom)
                    continue;
                end
            end
        end

        % Orientation
        orientation = ones(dofm.numDofPerTRI,1);
        if(ver(1) > ver(2))
            orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        if(ver(2) > ver(3))
            orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        if(ver(3) > ver(1))
            orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        orientation = sparse(1:dofm.numDofPerTRI, 1:dofm.numDofPerTRI, orientation);

        fieldTri = real(field(dofm.locToGloTRI(tri,:)));
        fieldTri = orientation*fieldTri;

        % Print
        fprintf(file,'%i %i ', mesh.tagTriFile(tri), dofm.numDofPerTRI);
        for n=1:dofm.numDofPerTRI
            fprintf(file,'%f ', fieldTri(n));
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
    fprintf(file,'%i\n',mesh.numTri-numTriPml);
    fprintf(file,'0\n');
    for tri=1:mesh.numTri
        ver = mesh.mapTriToVer(tri,:);

        if(PML_HIDE == 1)
            if(~isempty(LdomX) && ~isempty(LdomY))
                VX = mesh.coord(ver,1);
                VY = mesh.coord(ver,2);
                if ((mean(abs(VX)) >= LdomX) || (mean(abs(VY)) >= LdomY))
                    continue;
                end
            end
            if(~isempty(Rdom))
                VX = mesh.coord(ver,1);
                VY = mesh.coord(ver,2);
                VZ = VX + 1i*VY;
                VR = abs(VZ);
                if (mean(VR) >= Rdom)
                    continue;
                end
            end
        end

        % Orientation
        orientation = ones(dofm.numDofPerTRI,1);
        if(ver(1) > ver(2))
            orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        if(ver(2) > ver(3))
            orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        if(ver(3) > ver(1))
            orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        end
        orientation = sparse(1:dofm.numDofPerTRI, 1:dofm.numDofPerTRI, orientation);

        fieldTri = imag(field(dofm.locToGloTRI(tri,:)));
        fieldTri = orientation*fieldTri;

        % Print
        fprintf(file,'%i %i ', mesh.tagTriFile(tri), dofm.numDofPerTRI);
        for n=1:dofm.numDofPerTRI
            fprintf(file,'%f ', fieldTri(n));
        end
        fprintf(file,'\n');
    end
    fprintf(file,'$EndElementNodeData\n');

end

% Close file

fclose(file);

end
