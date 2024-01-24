% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function plotField1D(mesh, dofm, vecSol, nameTitle)

numVizuPntPerE = 30;
coordLoc  = linspace(-1,1,numVizuPntPerE)';
shapeFunc = functionsShapeLIN(coordLoc,dofm.degree);

vizuCoord = zeros(numVizuPntPerE,mesh.numE);
vizuVal   = zeros(numVizuPntPerE,mesh.numE);
for e=1:mesh.numE
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    vizuCoord(:,e) = linspace(coord1,coord2,numVizuPntPerE)';
    vizuVal(:,e) = shapeFunc * vecSol(dofm.locToGlo(e,:));
end

hold off
plot(vizuCoord(:),real(vizuVal(:)),'b','LineWidth',2,'DisplayName','Real part');
hold on
plot(vizuCoord(:),imag(vizuVal(:)),'r','LineWidth',2,'DisplayName','Imaginary part');
title(nameTitle);
lgd = legend;
lgd.FontSize = 12;

end