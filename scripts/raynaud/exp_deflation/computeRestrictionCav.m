function [matR] = computeRestrictionCav(mesh, dofm)


matR = sparse(dofm.numDofTRI,dofm.numDofTRI);

for tri=1:mesh.numTri

    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);

%     inX = V1(1) >= -0.75 && V2(1) >= -0.75 && V3(1) >= -0.75 && V1(1) <= 0.55  && V2(1) <= 0.55  && V3(1) <= 0.55;
%     inY = V1(2) >= -0.2  && V2(2) >= -0.2  && V3(2) >= -0.2  && V1(2) <= 0.2   && V2(2) <= 0.2   && V3(2) <= 0.2;

    inX = V1(1) > -0.75 && V2(1) > -0.75 && V3(1) > -0.75 && V1(1) <= 0.55  && V2(1) <= 0.55  && V3(1) <= 0.55;
    inY = V1(2) >= -0.2  && V2(2) >= -0.2  && V3(2) >= -0.2  && V1(2) <= 0.2   && V2(2) <= 0.2   && V3(2) <= 0.2;

    if inX && inY
        dof = dofm.locToGloTRI(tri,:);
        matR(dof,dof) = speye(numel(dof));
    end
end


end
