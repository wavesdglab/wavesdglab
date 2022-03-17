function postProVizu1D_DG(mesh, dofm, mySolP, mySolU, vecSol)

numVizuPntPerE = 50;
coordLoc = linspace(-1,1,numVizuPntPerE);
shapeFunc = functionsShape1D(coordLoc,dofm.degree);

vizuCoord  = [];
vizuValP = [];
vizuValU = [];
for e=1:mesh.numE
    
    coord1 = mesh.coordV(mesh.listE(e,1));
    coord2 = mesh.coordV(mesh.listE(e,2));
    length = abs(coord2-coord1);
    coordGlo = linspace(coord1,coord2,numVizuPntPerE);
    
    valP = zeros(1,numVizuPntPerE);
    valU = zeros(1,numVizuPntPerE);
    for n=1:dofm.numDofPerE
        valP = valP + vecSol(dofm.locToGlo(e,n)) * (shapeFunc(:,n))';
        valU = valU + vecSol(dofm.locToGlo(e,n) + dofm.numDof) * (shapeFunc(:,n))';
    end
    
    vizuCoord  = [vizuCoord coordGlo];
    vizuValP   = [vizuValP valP];
    vizuValU   = [vizuValU valU];
end

subplot(1,2,1);
hold off
plot(vizuCoord,real(vizuValP),'b','LineWidth',2,'DisplayName','Num solP (Real)');
title(['solP | numV = ' int2str(mesh.numV) ' | degree = ' int2str(dofm.degree)])
hold on
plot(vizuCoord,imag(vizuValP),'r','LineWidth',2,'DisplayName','Num solP (Imag)');
plot(vizuCoord,real(mySolP(vizuCoord)),':b','LineWidth',2.5,'DisplayName','Ref solP (Real)');
plot(vizuCoord,imag(mySolP(vizuCoord)),':r','LineWidth',2.5,'DisplayName','Ref solP (Real)');
%ylim([-1.5 1.5]);
%xlim([0 0.02]);

lgd = legend;
lgd.FontSize = 14;

subplot(1,2,2);
hold off
plot(vizuCoord,real(vizuValU),'b','LineWidth',2,'DisplayName','Num solU (Real)');
title(['solU | numV = ' int2str(mesh.numV) ' | degree = ' int2str(dofm.degree)])
hold on
plot(vizuCoord,imag(vizuValU),'r','LineWidth',2,'DisplayName','Num solU (Imag)');
plot(vizuCoord,real(mySolU(vizuCoord)),':b','LineWidth',2.5,'DisplayName','Ref solU (Real)');
plot(vizuCoord,imag(mySolU(vizuCoord)),':r','LineWidth',2.5,'DisplayName','Ref solU (Imag)');
%ylim([-1.5 1.5]);
%xlim([0 0.02]);
lgd = legend;
lgd.FontSize = 14;

end