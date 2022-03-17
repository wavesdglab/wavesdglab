function [matM, matK, matD] = buildMatrixGlo1D_CG(mesh, dofm)

[matElemM, matElemK, matElemD] = buildMatrixElem1D(dofm.degree);

matMgg = sparse(dofm.numDofGam, dofm.numDofGam);
matMgi = sparse(dofm.numDofGam, dofm.numDofInt);
matMig = sparse(dofm.numDofInt, dofm.numDofGam);
matMii = sparse(dofm.numDofInt, dofm.numDofInt);

matKgg = sparse(dofm.numDofGam, dofm.numDofGam);
matKgi = sparse(dofm.numDofGam, dofm.numDofInt);
matKig = sparse(dofm.numDofInt, dofm.numDofGam);
matKii = sparse(dofm.numDofInt, dofm.numDofInt);

matDgg = sparse(dofm.numDofGam, dofm.numDofGam);
matDgi = sparse(dofm.numDofGam, dofm.numDofInt);
matDig = sparse(dofm.numDofInt, dofm.numDofGam);
matDii = sparse(dofm.numDofInt, dofm.numDofInt);

for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2 - coord1);
    
    matLocM = matElemM * length/2;
    matLocK = matElemK * 2/length;
    matLocD = matElemD;
    
    glo = dofm.locToGlo(e,:);
    gloG = glo(1:2);
    gloI = glo(3:end)-dofm.numDofGam;
    
    matMgg(gloG,gloG) = matMgg(gloG,gloG) + matLocM(1:2,1:2);
    matMgi(gloG,gloI) = matMgi(gloG,gloI) + matLocM(1:2,3:end);
    matMig(gloI,gloG) = matMig(gloI,gloG) + matLocM(3:end,1:2);
    matMii(gloI,gloI) = matMii(gloI,gloI) + matLocM(3:end,3:end);
    
    matKgg(gloG,gloG) = matKgg(gloG,gloG) + matLocK(1:2,1:2);
    matKgi(gloG,gloI) = matKgi(gloG,gloI) + matLocK(1:2,3:end);
    matKig(gloI,gloG) = matKig(gloI,gloG) + matLocK(3:end,1:2);
    matKii(gloI,gloI) = matKii(gloI,gloI) + matLocK(3:end,3:end);
    
    matDgg(gloG,gloG) = matDgg(gloG,gloG) + matLocD(1:2,1:2);
    matDgi(gloG,gloI) = matDgi(gloG,gloI) + matLocD(1:2,3:end);
    matDig(gloI,gloG) = matDig(gloI,gloG) + matLocD(3:end,1:2);
    matDii(gloI,gloI) = matDii(gloI,gloI) + matLocD(3:end,3:end);
end

matM = [ matMgg matMgi ; matMig matMii ];
matK = [ matKgg matKgi ; matKig matKii ];
matD = [ matDgg matDgi ; matDig matDii ];

end