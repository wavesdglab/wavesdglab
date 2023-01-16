function [solPP, dofmPP] = computeSolPostPro2D_DG(mesh, dofm, solA)

global k

dofmPP = buildDofManager2D_DG(mesh, dofm.degree+1);

% Quadrature (degree+1)
degreeQ = 2*dofmPP.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);

% Shape functions
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);
shapePPQ = functionsShapeTRI(uQ, vQ, dofmPP.degree);
[shapePPDuQ, shapePPDvQ] = functionsShapeDerTRI(uQ, vQ, dofmPP.degree);

solP = solA((0*dofm.numDofTRI+1):(1*dofm.numDofTRI));
solU = solA((1*dofm.numDofTRI+1):(2*dofm.numDofTRI));
solV = solA((2*dofm.numDofTRI+1):(3*dofm.numDofTRI));
solPP = zeros(dofmPP.numDofTRI, 1);

for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    Jdudx = inv(Jdxdu);                 % [ du/dx du/dy ; dv/dx dv/dy ]
    detJdxdu = abs(det(Jdxdu));
    
    % Orientation
    orientation = ones(dofm.numDofPerTRI,1);
    orientationPP = ones(dofmPP.numDofPerTRI,1);
    if(ver(1) > ver(2))
        orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        orientationPP(dofmPP.locEdg(1,:)) = (-1).^(0:dofmPP.numDofPerEdg-1);
    end
    if(ver(2) > ver(3))
        orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        orientationPP(dofmPP.locEdg(2,:)) = (-1).^(0:dofmPP.numDofPerEdg-1);
    end
    if(ver(3) > ver(1))
        orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
        orientationPP(dofmPP.locEdg(3,:)) = (-1).^(0:dofmPP.numDofPerEdg-1);
    end
    orientation = sparse(1:dofm.numDofPerTRI, 1:dofm.numDofPerTRI, orientation);
    orientationPP = sparse(1:dofmPP.numDofPerTRI, 1:dofmPP.numDofPerTRI, orientationPP);
    
    % Shape functions with orientation
    shapeOrQ = shapeQ * orientation;
    shapePPOrQ = shapePPQ * orientationPP;
    shapePPDxQ = (shapePPDuQ * Jdudx(1,1) + shapePPDvQ * Jdudx(2,1)) * orientationPP;
    shapePPDyQ = (shapePPDuQ * Jdudx(1,2) + shapePPDvQ * Jdudx(2,2)) * orientationPP;
    
    % DOF
    dof = dofm.locToGloTRI(tri,:);
    solPel = solP(dof);
    solUel = solU(dof);
    solVel = solV(dof);
    
    % Elemental system
    matKel = shapePPDxQ' * weights * shapePPDxQ * detJdxdu + shapePPDyQ' * weights * shapePPDyQ * detJdxdu;
    rhsKel = 1i*k * (shapePPDxQ' * weights * shapeOrQ * solUel + shapePPDyQ' * weights * shapeOrQ * solVel) * detJdxdu;
    matCel = ones(1,size(shapePPOrQ,1)) * weights * shapePPOrQ * detJdxdu;
    rhsCel = ones(1,size(shapeOrQ,1)) * weights * shapeOrQ * solPel * detJdxdu;
    
    matPPel = [ matKel matCel' ; matCel 0 ];
    rhsPPel = [ rhsKel ; rhsCel ];
    solPPel = matPPel \ rhsPPel;
    
    % Matrix assembling
    dofPP = dofmPP.locToGloTRI(tri,:);
    solPP(dofPP) = solPPel(1:dofmPP.numDofPerTRI);
    
end

end