function dofm = buildDofManager2D_DG(mesh, degree)

% Polynomial degree
dofm.degree = degree;

% Number of DOF per vertex, edge and face
dofm.numDofPerVer = 1;
dofm.numDofPerEdg = dofm.degree-1;
dofm.numDofPerFac = max(0, (dofm.degree-1)*(dofm.degree-2)/2);

% Number of DOF per LIN and TRI element
dofm.numDofPerLIN = 2*dofm.numDofPerVer + dofm.numDofPerEdg;
dofm.numDofPerTRI = 3*dofm.numDofPerVer + 3*dofm.numDofPerEdg + dofm.numDofPerFac;

% Total number of DOF for LIN and TRI elements on the mesh
dofm.numDofLIN = dofm.numDofPerLIN * mesh.numEdg;
dofm.numDofTRI = dofm.numDofPerTRI * mesh.numTri;

% Mapping element-local to mesh-global index of DOF
dofm.locToGloLIN = zeros(mesh.numEdg, dofm.numDofPerLIN);
dofm.locToGloTRI = zeros(mesh.numTri, dofm.numDofPerTRI);
for edg=1:mesh.numEdg
    for n=1:dofm.numDofPerLIN
        dofm.locToGloLIN(edg,n) = (edg-1)*dofm.numDofPerLIN + n;
    end
end
for tri=1:mesh.numTri
    for n=1:dofm.numDofPerTRI
        dofm.locToGloTRI(tri,n) = (tri-1)*dofm.numDofPerTRI + n;
    end
end

end