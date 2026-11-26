% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function plotMesh2D(mesh)

figure;
trimesh(mesh.mapTriToVer, mesh.coord(:,1), mesh.coord(:,2), zeros(mesh.numVer,1));
view(2);
axis('equal');

end