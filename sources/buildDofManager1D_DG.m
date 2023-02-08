% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function dofm = buildDofManager1D_DG(mesh, degree)

% Polynomial degree
dofm.degree = degree;

% Number of DOF per element
dofm.numDofPerE = degree+1;

% Total number of DOF on the mesh
dofm.numDof = mesh.numE * dofm.numDofPerE;

% Mapping element-local to mesh-global index of DOF
dofm.locToGlo = zeros(mesh.numE, dofm.numDofPerE);
for e=1:mesh.numE
    for d=1:dofm.numDofPerE
        dofm.locToGlo(e,d) = (e-1)*dofm.numDofPerE + d;
    end
end

end