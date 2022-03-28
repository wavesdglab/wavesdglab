function [normL2sol, normL2der, normH1sol] = computeNormSol2D(mesh, degreeQ)

[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);

normL2sol2 = 0;
normL2der2 = 0;
for tri=1:mesh.numTri
    
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uQ, vQ, V1, V2, V3);
    
    [solQ, solDxQ, solDyQ] = mySol(xQ, yQ);
    
    Jdxdu = [(V3-V2)' (V1-V2)'] * 0.5;
    detJdxdu = abs(det(Jdxdu));
    
    normL2sol2 = normL2sol2 + weights(:)' * (solQ .* conj(solQ)) * detJdxdu;
    normL2der2 = normL2der2 + weights(:)' * (solDxQ .* conj(solDxQ) + solDyQ .* conj(solDyQ)) * detJdxdu;
end

normL2sol = sqrt(normL2sol2);
normL2der = sqrt(normL2der2);
normH1sol = sqrt(normL2sol2 + normL2der2);

end