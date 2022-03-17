function postProVizu1D_CG(mesh, dofm, mySol, mySolDer, vecSol)

numVizuPntPerE = 30;
coordLoc = linspace(-1,1,numVizuPntPerE);
shapeFunc    = functionsShape1D(coordLoc,dofm.degree);
shapeFuncDer = functionsShapeDer1D(coordLoc,dofm.degree);

% figure(3);
% plot(shapeFunc,'LineWidth',2);
% 
% figure(4);
% plot(shapeFuncDer,'LineWidth',2);

vizuCoord  = [];
vizuVal    = [];
vizuValDer = [];
for e=1:mesh.numE
    
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2-coord1);
    coordGlo = linspace(coord1,coord2,numVizuPntPerE);
    
    val = zeros(1,numVizuPntPerE);
    valDer = zeros(1,numVizuPntPerE);
    for n=1:dofm.numDofPerE
        val    = val    + vecSol(dofm.locToGlo(e,n)) * (shapeFunc(:,n)   )';
        valDer = valDer + vecSol(dofm.locToGlo(e,n)) * (shapeFuncDer(:,n))' * (2/length);
    end
    
    vizuCoord  = [vizuCoord coordGlo];
    vizuVal    = [vizuVal val];
    vizuValDer = [vizuValDer valDer];
end

hold off
plot(vizuCoord,real(vizuVal),'b','LineWidth',2,'DisplayName','Num sol (Real)');
title(['sol | numV = ' int2str(mesh.numV) ' | degree = ' int2str(dofm.degree)])
hold on
plot(vizuCoord,imag(vizuVal),'r','LineWidth',2,'DisplayName','Num sol (Imag)');
plot(vizuCoord,real(mySol(vizuCoord)),':b','LineWidth',2.5,'DisplayName','Ref sol (Real)');
plot(vizuCoord,imag(mySol(vizuCoord)),':r','LineWidth',2.5,'DisplayName','Ref sol (Real)');

lgd = legend;
lgd.FontSize = 14;

% subplot(1,2,2);
% % hold off
% % plot(vizuCoord,real(vizuVal)-real(mySol(vizuCoord)));
% % title(['error | numV = ' int2str(numV) ' | degree = ' int2str(degree)])
% % hold on
% % plot(vizuCoord,imag(vizuVal)-imag(mySol(vizuCoord)));
% hold off
% plot(vizuCoord,real(vizuValDer),'b','LineWidth',2,'DisplayName','Num derSol (Real)');
% title(['derSol | numV = ' int2str(mesh.numV) ' | degree = ' int2str(dofm.degree)])
% hold on
% plot(vizuCoord,imag(vizuValDer),'r','LineWidth',2,'DisplayName','Num derSol (Imag)');
% plot(vizuCoord,real(mySolDer(vizuCoord)),':b','LineWidth',2.5,'DisplayName','Ref derSol (Real)');
% plot(vizuCoord,imag(mySolDer(vizuCoord)),':r','LineWidth',2.5,'DisplayName','Ref derSol (Imag)');
% 
% lgd = legend;
% lgd.FontSize = 14;

end