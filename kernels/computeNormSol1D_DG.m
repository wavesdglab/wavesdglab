function normL2sol = computeNormSol1D_DG(mesh, mySol)

Q = 16;
[nodes, weights] = quadratureGauss1D(Q);

normL2sol2 = 0;
for e=1:mesh.numE
    
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2-coord1);
    coordGlo = coord1*(1-nodes)/2 + coord2*(1+nodes)/2;
    
    valRef = mySol(coordGlo);
    
    normL2sol2 = normL2sol2 + weights(:)' * (valRef.*conj(valRef)) * (length/2) ;
end

normL2sol = sqrt(normL2sol2);

end