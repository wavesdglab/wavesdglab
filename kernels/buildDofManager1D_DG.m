function dofm = buildDofManager1D_DG(mesh, degree)

dofm.degree     = degree;
dofm.numDofPerE = degree+1;
dofm.numDof     = mesh.numE * dofm.numDofPerE;

dofm.locToGlo = zeros(mesh.numE, dofm.numDofPerE);
for e=1:mesh.numE
    for d=1:dofm.numDofPerE
        dofm.locToGlo(e,d) = (e-1)*dofm.numDofPerE + d;
    end
end

end