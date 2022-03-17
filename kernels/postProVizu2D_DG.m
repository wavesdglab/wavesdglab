function postProVizu2D_DG(mesh, field, titre)

for tri=1:mesh.numTri
    nodGloTri = 3*(tri-1)+1:3*tri;
    trisurf(1:3, mesh.coord(mesh.mapTriToVer(tri,:),1), mesh.coord(mesh.mapTriToVer(tri,:),2), field(nodGloTri));
    hold on;
end
view(2);
shading interp
colorbar;
title(titre);
