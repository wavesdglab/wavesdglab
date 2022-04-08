function dofm = buildDofManager2D_CG(mesh, degree)

% Polynomial degree
dofm.degree = degree;

% Number of DOF per vertex, edge and face
dofm.numDofPerVer = 1;
dofm.numDofPerEdg = dofm.degree-1;
dofm.numDofPerFac = (dofm.degree-1)*(dofm.degree-2)/2;

% Number of DOF per LIN and TRI element
dofm.numDofPerLIN = 2*dofm.numDofPerVer + dofm.numDofPerEdg;
dofm.numDofPerTRI = 3*dofm.numDofPerVer + 3*dofm.numDofPerEdg + dofm.numDofPerFac;

% Total number of DOF associated to vertices, edges and faces on the mesh
dofm.numDofVer = mesh.numVer * dofm.numDofPerVer;
dofm.numDofEdg = mesh.numEdg * dofm.numDofPerEdg;
dofm.numDofFac = mesh.numTri * dofm.numDofPerFac;

% Total number of DOF on the mesh
dofm.numDofTRI = dofm.numDofVer + dofm.numDofEdg + dofm.numDofFac;

% Mapping element-local to mesh-global index of vertex/edge/face DOF
dofm.locToGloVer = zeros(mesh.numTri, 3*dofm.numDofPerVer);
dofm.locToGloEdg = zeros(mesh.numTri, 3*dofm.numDofPerEdg);
dofm.locToGloFac = zeros(mesh.numTri, 1*dofm.numDofPerFac);
dofm.locEdg = zeros(3,dofm.numDofPerEdg);
for tri=1:mesh.numTri
    ver = mesh.mapTriToVer(tri,:);
    dofm.locToGloVer(tri,1) = ver(1);
    dofm.locToGloVer(tri,2) = ver(2);
    dofm.locToGloVer(tri,3) = ver(3);
    edg = abs(mesh.mapTriToEdg(tri,:));
    for n=1:dofm.numDofPerEdg
        dofm.locToGloEdg(tri,0*dofm.numDofPerEdg+n) = dofm.numDofVer + dofm.numDofPerEdg*(edg(1)-1) + n;
        dofm.locToGloEdg(tri,1*dofm.numDofPerEdg+n) = dofm.numDofVer + dofm.numDofPerEdg*(edg(2)-1) + n;
        dofm.locToGloEdg(tri,2*dofm.numDofPerEdg+n) = dofm.numDofVer + dofm.numDofPerEdg*(edg(3)-1) + n;
        dofm.locEdg(1,n) = 3 + 0*dofm.numDofPerEdg + n;
        dofm.locEdg(2,n) = 3 + 1*dofm.numDofPerEdg + n;
        dofm.locEdg(3,n) = 3 + 2*dofm.numDofPerEdg + n;
    end
    for n=1:dofm.numDofPerFac
        dofm.locToGloFac(tri,n) = dofm.numDofVer + dofm.numDofEdg + dofm.numDofPerFac*(tri-1) + n;
    end
end
dofm.locToGloTRI = [dofm.locToGloVer dofm.locToGloEdg dofm.locToGloFac];

% Mapping element-local to mesh-global index of vertex/edge/face DOF
dofm.locBndToGloVer = zeros(mesh.numEdgBnd, 2*dofm.numDofPerVer);
dofm.locBndToGloEdg = zeros(mesh.numEdgBnd, dofm.numDofPerEdg);
for edgBnd=1:mesh.numEdgBnd
    edg = mesh.listEdgBnd(edgBnd);
    ver = mesh.mapEdgToVer(edg,:);
    dofm.locBndToGloVer(edgBnd,1) = ver(1);
    dofm.locBndToGloVer(edgBnd,2) = ver(2);
    for n=1:dofm.numDofPerEdg
        dofm.locBndToGloEdg(edgBnd,n) = dofm.numDofVer + dofm.numDofPerEdg*(edg-1) + n;
    end
end
dofm.locToGloBND = [dofm.locBndToGloVer dofm.locBndToGloEdg];

end