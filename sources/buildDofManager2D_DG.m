% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function dofm = buildDofManager2D_DG(mesh, degree)

% Polynomial degree
dofm.degree = degree;

% Number of DOF per vertex, edge and face
dofm.numDofPerVer = 1;
dofm.numDofPerEdg = dofm.degree-1;
dofm.numDofPerFac = (dofm.degree-1)*(dofm.degree-2)/2;

% Number of DOF per LIN and TRI element
dofm.numDofPerLIN = 2*dofm.numDofPerVer + dofm.numDofPerEdg;
dofm.numDofPerTRI = 3*dofm.numDofPerVer + 3*dofm.numDofPerEdg + dofm.numDofPerFac;

% Total number of DOF for LIN and TRI elements on the mesh
dofm.numDofLIN = dofm.numDofPerLIN * mesh.numEdg;
dofm.numDofTRI = dofm.numDofPerTRI * mesh.numTri;
dofm.numDofFAC = dofm.numDofPerLIN * mesh.numTri * 3;

% Mapping edge-local to element-local index of edge/face DOF
dofm.locEdg = zeros(3,dofm.numDofPerEdg);
dofm.locEdg(1,:) = 3 + 0*dofm.numDofPerEdg + (1:dofm.numDofPerEdg);
dofm.locEdg(2,:) = 3 + 1*dofm.numDofPerEdg + (1:dofm.numDofPerEdg);
dofm.locEdg(3,:) = 3 + 2*dofm.numDofPerEdg + (1:dofm.numDofPerEdg);
dofm.locFac = zeros(3,dofm.numDofPerLIN);
dofm.locFac(1,:) = [1 2 dofm.locEdg(1,:)];
dofm.locFac(2,:) = [2 3 dofm.locEdg(2,:)];
dofm.locFac(3,:) = [3 1 dofm.locEdg(3,:)];
dofm.locFacNeigh = zeros(3,dofm.numDofPerLIN);
dofm.locFacNeigh(1,:) = [2 1 dofm.locEdg(1,:)];
dofm.locFacNeigh(2,:) = [3 2 dofm.locEdg(2,:)];
dofm.locFacNeigh(3,:) = [1 3 dofm.locEdg(3,:)];

% Mapping element-local to mesh-global index of LIN/TRI DOF
dofm.locToGloLIN = zeros(mesh.numEdg, dofm.numDofPerLIN);
dofm.locToGloTRI = zeros(mesh.numTri, dofm.numDofPerTRI);
dofm.locToGloFAC = zeros(mesh.numTri, 3*dofm.numDofPerLIN);
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
for tri=1:mesh.numTri
    for n=1:3*dofm.numDofPerLIN
        dofm.locToGloFAC(tri,n) = (tri-1)*3*dofm.numDofPerLIN + n;
    end
end

end