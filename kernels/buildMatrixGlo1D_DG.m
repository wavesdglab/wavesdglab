function [matM, matK, matD] = buildMatrixGlo1D_DG(mesh, dofm)

[matElemM, matElemK, matElemD] = buildMatrixElem1D(dofm.degree);

matM = sparse(dofm.numDof, dofm.numDof);
matK = sparse(dofm.numDof, dofm.numDof);
matD = sparse(dofm.numDof, dofm.numDof);

for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    
    matLocM = matElemM * length/2;
    matLocK = matElemK * 2/length;
    matLocD = matElemD;
    
    glo = dofm.locToGlo(e,:);
    matM(glo,glo) = matM(glo,glo) + matLocM;
    matK(glo,glo) = matK(glo,glo) + matLocK;
    matD(glo,glo) = matD(glo,glo) + matLocD;
end

end