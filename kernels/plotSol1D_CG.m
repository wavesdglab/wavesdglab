function plotSol1D_CG(mesh, dofm, vecSol, nameTitle)

numVizuPntPerE = 30;
coordLoc  = linspace(-1,1,numVizuPntPerE);
shapeFunc = functionsShape1D(coordLoc,dofm.degree);

vizuCoord = [];
vizuVal   = [];
for e=1:mesh.numE
    
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    coordGlo = linspace(coord1,coord2,numVizuPntPerE);
    
    val = zeros(1,numVizuPntPerE);
    for n=1:dofm.numDofPerE
        val = val + vecSol(dofm.locToGlo(e,n)) * (shapeFunc(:,n))';
    end
    
    vizuCoord = [vizuCoord coordGlo];
    vizuVal   = [vizuVal val];
end

hold off
plot(vizuCoord,real(vizuVal),'b','LineWidth',2,'DisplayName','Real part');
hold on
plot(vizuCoord,imag(vizuVal),'r','LineWidth',2,'DisplayName','Imaginary part');
title(nameTitle);
lgd = legend;
lgd.FontSize = 14;

end