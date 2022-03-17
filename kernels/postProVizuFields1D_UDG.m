function postProVizuFields1D_UDG(mesh, dofm, vecSol)

numVizuPntPerE = 30;
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
plot(vizuCoord,real(vizuValP),'LineWidth',2);
subplot(1,2,2);
plot(vizuCoord,imag(vizuValP),'LineWidth',2);

% subplot(1,2,1);
% hold off
% plot(vizuCoord,real(vizuValP),'b','LineWidth',2);
% hold on
% plot(vizuCoord,imag(vizuValP),'r','LineWidth',2);
% xlim([0 1]);
% ylim([-1 1]);
% 
% subplot(1,2,2);
% hold off
% plot(vizuCoord,real(vizuValU),'b','LineWidth',2);
% hold on
% plot(vizuCoord,imag(vizuValU),'r','LineWidth',2);
% xlim([0 1]);
% ylim([-1 1]);

end