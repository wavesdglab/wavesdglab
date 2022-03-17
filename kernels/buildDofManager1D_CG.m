function dofm = buildDofManager1D_CG(mesh, degree)

global BCLeft BCRight

dofm.degree        = degree;
dofm.numDofPerE    = degree+1;
dofm.numDofPerEInt = degree-1;

dofm.numDof    = mesh.numV + mesh.numE * dofm.numDofPerEInt;
dofm.numDofGam = mesh.numV;
dofm.numDofInt = mesh.numE * dofm.numDofPerEInt;

if(strcmp(BCLeft,'PER') && strcmp(BCRight,'PER'))
    dofm.numDof = dofm.numDof-1;
    dofm.numDofGam = dofm.numDofGam-1;
end

dofm.locToGlo = zeros(mesh.numE, dofm.numDofPerE);
for e=1:mesh.numE
    dofm.locToGlo(e,1) = e;
    dofm.locToGlo(e,2) = e+1;
    for d=3:dofm.numDofPerE
        dofm.locToGlo(e,d) = dofm.numDofGam + (e-1)*dofm.numDofPerEInt + (d-2);
    end
end

if(strcmp(BCLeft,'PER') && strcmp(BCRight,'PER'))
    dofm.locToGlo(mesh.numE,2) = dofm.locToGlo(1,1);
end

end