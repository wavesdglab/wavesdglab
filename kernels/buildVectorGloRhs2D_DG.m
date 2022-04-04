function [rhsGlo] = buildVectorGloRhs2D_DG(mesh, dofm, rhsFun)

% Quadrature
degreeQ = 4*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);

% Shape functions
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);

% Global matrices
rhsGlo = zeros(dofm.numDofTRI,1);

for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uQ, vQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    detJdxdu = abs(det(Jdxdu));
    
    % RHS function
    rhsQ = rhsFun(xQ, yQ);
    
    % Elemental matrices
    rhsLoc = zeros(dofm.numDofPerTRI,1);
    for i=1:dofm.numDofPerTRI
        rhsLoc(i) = weights' * (shapeQ(:,i) .* rhsQ(:)) * detJdxdu;
    end
    
    % Matrix assembling
    dof = dofm.locToGloTRI(tri,:);
    rhsGlo(dof) = rhsLoc;
    
end

end