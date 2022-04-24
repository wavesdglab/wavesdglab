function [normL2sol, normL2der, normH1sol] = computeNormSol1D_CG(mesh, mySol, mySolDer)

Q = 16;
[nodes, weights] = quadratureGaussLIN(Q);

normL2sol2 = 0;
normL2der2 = 0;
for e=1:mesh.numE
    
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2-coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    solQ = mySol(coordGlo);
    solDerQ = mySolDer(coordGlo);
    
    normL2sol2 = normL2sol2 + weights(:)' * (solQ .* conj(solQ))       * (length/2);
    normL2der2 = normL2der2 + weights(:)' * (solDerQ .* conj(solDerQ)) * (length/2);
end

normL2sol = sqrt(normL2sol2);
normL2der = sqrt(normL2der2);
normH1sol = sqrt(normL2sol2 + normL2der2);

end