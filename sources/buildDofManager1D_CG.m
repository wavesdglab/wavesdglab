function dofm = buildDofManager1D_CG(mesh, degree)

global BCLeft BCRight

% Polynomial degree
dofm.degree = degree;

% Number of DOF per element
dofm.numDofPerE    = degree+1;  % All DOF
dofm.numDofPerEInt = degree-1;  % Only interior DOF

% Total number of DOF associated to vertices and edges on the mesh
dofm.numDofGam = mesh.numV;                       % associated to vertices
dofm.numDofInt = mesh.numE * dofm.numDofPerEInt;  % associated to edges

if(strcmp(BCLeft,'PER') && strcmp(BCRight,'PER'))
    dofm.numDofGam = dofm.numDofGam-1;
end

% Total number of DOF on the mesh
dofm.numDof = dofm.numDofGam + dofm.numDofInt;

% Mapping element-local to mesh-global index of DOF
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