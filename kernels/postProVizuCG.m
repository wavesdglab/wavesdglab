function postProVizuCG(mesh, field, titre)

trisurf(mesh.mapTriToVer, mesh.coord(:,1), mesh.coord(:,2), field);
view(2);
shading interp
colorbar;
title(titre);
