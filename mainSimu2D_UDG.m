%close all;
%clear all;

headers2D;

% Define parameters
global k BCWest BCNorth BCEast BCSouth
k = 3;
h = 0.1;
degree = 3;
tau = 1;
resTol = 1e-4;
BCWest  = 'ABC';
BCNorth = 'DIR';
BCEast  = 'NEU';
BCSouth = 'ABC';

% Build mesh and dofManager
system(['gmsh -2 mesh.geo -v 0 -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh('mesh.msh');
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[solA, matA, rhsA] = computeSolNum2D_UDG1(mesh, dofm, tau);
[errorL2, errorH1] = computeNormError2D_DG(mesh, dofm, solA);

[solP, matP, rhsP] = computeSolProjL2_2D_DG(mesh, dofm);
[errorProjL2, errorProjH1] = computeNormError2D_DG(mesh, dofm, solP);

disp('Method HDG-1');
disp('---------------------------------------------------------');
disp(['    L2-Error (numSol)  ' num2str(errorL2)]);
disp(['    H1-Error (numSol)  ' num2str(errorH1)]);
disp(['    L2-Error (projSol) ' num2str(errorProjL2)]);
disp(['    H1-Error (projSol) ' num2str(errorProjH1)]);
disp('---------------------------------------------------------');

%writeFieldDG(dofm, mesh, solP, "mySol.pos", "mySol");
%system('gmsh mySol.pos');


% fprintf('Solver  : gmres A\n');
% [solAiter,~,~,iterA]         = gmres(matA,rhsA,size(matA,1),resTol,size(matA,1));
% errorL2IterA             = computeNormError2D_DG(mesh, dofm, solAiter, solA);
% fprintf('Solver  : bicgstab A\n');
% [solAiter,~,~,iterBiCGStabA] = bicgstab(matA,rhsA,resTol,size(matA,1));
% errorL2BiCGStabA         = computeNormError2D_DG(mesh, dofm, solAiter, solA);
% fprintf('Solver  : conjgradn A\n');
% [solAiter,~,~,iterCGNA]      = conjgradn(matA,rhsA,resTol,size(matA,1));
% errorL2CGNA              = computeNormError2D_DG(mesh, dofm, solAiter, solA);
% fprintf('Solver  : relaxation A\n');
% [solAiter,~,~,iterJacobiA]   = jacobi(matA,rhsA,resTol,size(matA,1),0.5);
% errorL2JacobiA           = computeNormError2D_DG(mesh, dofm, solAiter, solA);

% matA11 = matA(1:3*dofm.numDofTRI, 1:3*dofm.numDofTRI);
% matA12 = matA(1:3*dofm.numDofTRI, 3*dofm.numDofTRI+1:end);
% matA21 = matA(3*dofm.numDofTRI+1:end, 1:3*dofm.numDofTRI);
% matA22 = matA(3*dofm.numDofTRI+1:end, 3*dofm.numDofTRI+1:end);
% rhsA1 = rhsA(1:3*dofm.numDofTRI);
% rhsA2 = rhsA(3*dofm.numDofTRI+1:end);
% matS = matA22 - matA21*(matA11\matA12);
% rhsS = rhsA2 - matA21*(matA11\rhsA1);
% 
% fprintf('Solver  : gmres S\n');
% [solS,~,~,iterS]         = gmres(matS,rhsS,size(matS,1),resTol,size(matS,1));
% solAiter = matA11\(rhsA1 - matA12*solS);
% errorL2IterS             = computeNormError2D_DG(mesh, dofm, solAiter, solA);
% fprintf('Solver  : bicgstab S\n');
% [solS,~,~,iterBiCGStabS] = bicgstab(matS,rhsS,resTol,size(matS,1));
% solAiter = matA11\(rhsA1 - matA12*solS);
% errorL2BiCGStabS         = computeNormError2D_DG(mesh, dofm, solAiter, solA);
% fprintf('Solver  : conjgradn S\n');
% [solS,flag,~,iterCGNS]      = conjgradn(matS,rhsS,resTol,size(matS,1));
% solAiter = matA11\(rhsA1 - matA12*solS);
% errorL2CGNS              = computeNormError2D_DG(mesh, dofm, solAiter, solA);
% [solS,flag,~,iterCGNS]      = conjgradnResPhy(matS,rhsS,resTol,size(matS,1),matA11,matA12,matA21,matA22,rhsA1,rhsA2);
% solAiter = matA11\(rhsA1 - matA12*solS);
% errorL2CGNS              = computeNormError2D_DG(mesh, dofm, solAiter, solA);
% fprintf('Solver  : relaxation S\n');
% [solS,~,~,iterJacobiS]   = jacobi(matS,rhsS,resTol,size(matS,1),0.5);
% solAiter = matA11\(rhsA1 - matA12*solS);
% errorL2JacobiS           = computeNormError2D_DG(mesh, dofm, solAiter, solA);
% fprintf('Solver  : richardson S\n');
% [solS,~,~,iterRichS]   = richardsonResPhy(matS,rhsS,resTol,size(matS,1),0.5,matA11,matA12,matA21,matA22,rhsA1,rhsA2);
% solAiter = matA11\(rhsA1 - matA12*solS);
% errorL2RichS           = computeNormError2D_DG(mesh, dofm, solAiter, solA);

% fprintf('Solver  : eigenval A\n');
% [eigenvecA,eigenvalA] = eigs(matA,size(matA,1));
% eigenvalA = diag(eigenvalA);
% [eigenvecAA,eigenvalAA] = eigs(matA'*matA,size(matA,1));
% eigenvalAA = diag(eigenvalAA);
% 
% fprintf('Solver  : eigenval S\n');
% [eigenvecS,eigenvalS] = eigs(matS,size(matS,1));
% eigenvalS = diag(eigenvalS);
% [eigenvecSS,eigenvalSS] = eigs(matS'*matS,size(matS,1));
% eigenvalSS = diag(eigenvalSS);

% disp(['A : Size               ' num2str(size(matA,1))]);
% disp(['    Rank(eigenvectors) ' num2str(rank(eigenvecA))]);
% disp(['    Cond(eigenvectors) ' num2str(cond(eigenvecA))]);
% disp(['    Cond(A)            ' num2str(condest(matA))]);
% disp(['    Cond(AA)           ' num2str(condest(matA'*matA))]);
% disp(['    Min eigenval AA    ' num2str(min(eigenvalAA))]);
% disp(['    - - - - - - - - - - - - - - - - - - - - -']);
% disp(['    IterGmres          ' num2str(iterA(2))]);
% disp(['    Final L2-Error     ' num2str(errorL2IterA)]);
% disp(['    IterBiCGS          ' num2str(iterBiCGStabA)]);
% disp(['    Final L2-Error     ' num2str(errorL2BiCGStabA)]);
% disp(['    IterCGN            ' num2str(iterCGNA)]);
% disp(['    Final L2-Error     ' num2str(errorL2CGNA)]);
% disp(['    IterRelax          ' num2str(iterJacobiA)]);
% disp(['    Final L2-Error     ' num2str(errorL2JacobiA)]);
% disp(['---------------------------------------------------------']);
% disp(['S : Size               ' num2str(size(matS,1))]);
% disp(['    Rank(eigenvectors) ' num2str(rank(eigenvecS))]);
% disp(['    Cond(eigenvectors) ' num2str(cond(eigenvecS))]);
% disp(['    Cond(S)            ' num2str(condest(matS))]);
% disp(['    Cond(SS)           ' num2str(condest(matS'*matS))]);
% disp(['    Min eigenval SS    ' num2str(min(eigenvalS))]);
% disp(['    - - - - - - - - - - - - - - - - - - - - -']);
% disp(['    IterGmres          ' num2str(iterS(2))]);
% disp(['    Final L2-Error     ' num2str(errorL2IterS)]);
% disp(['    IterBiCGS          ' num2str(iterBiCGStabS)]);
% disp(['    Final L2-Error     ' num2str(errorL2BiCGStabS)]);
% disp(['    IterCGN            ' num2str(iterCGNS)]);
% disp(['    Final L2-Error     ' num2str(errorL2CGNS)]);
% disp(['    IterRelax          ' num2str(iterJacobiS)]);
% disp(['    Final L2-Error     ' num2str(errorL2JacobiS)]);
% disp(['    IterRichard        ' num2str(iterRichS)]);
% disp(['    Final L2-Error     ' num2str(errorL2RichS)]);
% disp(['---------------------------------------------------------']);

% disp(['\text{UDG-1}(\tau=1) & full & ' ...
%     num2str(errorL2,'%.1e') ' & ' ...
%     num2str(errorH1,'%.1e') ' & ' ...
%     num2str(size(matA,1)) ' & ' ...
%     num2str(rank(eigenvecA)) ' & ' ...
%     num2str(cond(eigenvecA),'%.1e') ' & ' ...
%     num2str(condest(matA),'%.1e') ' & ' ...
%     num2str(iterA(2)) ' & ' ...
%     num2str(errorL2IterA,'%.1e') ' & ' ...
%     num2str(iterBiCGStabA) ' & ' ...
%     num2str(errorL2BiCGStabA,'%.1e') ' & ' ...
%     num2str(iterCGNA) ' & ' ...
%     num2str(errorL2CGNA,'%.1e') ' & ' ...
%     num2str(iterJacobiA) ' & ' ...
%     num2str(errorL2JacobiA,'%.1e') ' \\'
%     ]);

% disp(['\text{UDG-1}(\tau=i) & red & ' ...
%     num2str(errorL2,'%.1e') ' & ' ...
%     num2str(errorH1,'%.1e') ' & ' ...
%     num2str(size(matS,1)) ' & ' ...
%     num2str(rank(eigenvecS)) ' & ' ...
%     num2str(cond(eigenvecS),'%.1e') ' & ' ...
%     num2str(condest(matS),'%.1e') ' & ' ...
%     num2str(iterS(2)) ' & ' ...
%     num2str(errorL2IterS,'%.1e') ' & ' ...
%     num2str(iterBiCGStabS) ' & ' ...
%     num2str(errorL2BiCGStabS,'%.1e') ' & ' ...
%     num2str(iterCGNS) ' & ' ...
%     num2str(errorL2CGNS,'%.1e') ' & ' ...
%     num2str(iterJacobiS) ' & ' ...
%     num2str(errorL2JacobiS,'%.1e') ' \\'
%     ]);

% figure;
% hold off
% scatter(real(eigenvalS),imag(eigenvalS),'b','DisplayName','Eigenvalues');
% hold on
% %plot(fovals(matS,100),'-b','DisplayName','Numerical range');
% grid on; box on;
% title(['Eigenvalues : ' BCWest ' + ' BCNorth ' + ' BCEast ' + ' BCSouth]);
% legend();
% %axis([-0.1 1.1 -1.5 1.2]);
% %axis([-0.05 0.5 -1.5 1.2]);
% 
% figure(1);
% subplot(1,3,1);
% hold off
% postProVizu2D_DG(mesh,real(solP), 'Exact solution');
% subplot(1,3,2);
% hold off
% postProVizu2D_DG(mesh,real(solA), 'Numerical solution');
% subplot(1,3,3);
% hold off
% postProVizu2D_DG(mesh,real(solA-solP), 'Error');