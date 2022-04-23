close all;
clear all;

headers2D;

% Define parameters
global k
k = 20;
h = 0.1;
degree = 3;
tau = 1i;
resTol = 1e-4;
benchmark2D('cavity');

% Build mesh and dofManager
system(['gmsh -2 mesh.geo -v 0 -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh('mesh.msh');
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method HDG-1']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_HDG1(mesh, dofm, tau);
[errorL2, errorH1] = computeNormError2D_DG(mesh, dofm, solA);

[solP, matP, rhsP] = computeSolProjL2_2D_DG(mesh, dofm);
[errorProjL2, errorProjH1] = computeNormError2D_DG(mesh, dofm, solP);

[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solA);
[errorPostL2, errorPostH1] = computeNormError2D_DG(mesh, dofmPost, solApost);

[solPpost, matPpost, rhsPpost] = computeSolProjL2_2D_DG(mesh, dofmPost);
[errorProjPostL2, errorProjPostH1] = computeNormError2D_DG(mesh, dofmPost, solPpost);

disp(['---------------------------------------------------------']);
disp(['    L2-Error (numSol)   ' num2str(errorL2)]);
%disp(['    H1-Error (numSol)   ' num2str(errorH1)]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2)]);
%disp(['    H1-Error (projSol)  ' num2str(errorProjH1)]);
disp(['    L2-Error (numPost)  ' num2str(errorPostL2)]);
%disp(['    H1-Error (numPost)  ' num2str(errorPostH1)]);
disp(['    L2-Error (projPost) ' num2str(errorProjPostL2)]);
%disp(['    H1-Error (projPost) ' num2str(errorProjPostH1)]);
disp('---------------------------------------------------------');

% ValH = [0.2 0.1 0.05];
% ValSans = [0.04475 0.0020197 0.00023045];   % 8.76
% ValAvec = [0.039575 0.0010212 0.00011511];  % 8.87
% ValCea  = [0.019227 0.0015966 0.00020613];  % 7.74
% loglog(1./ValH, ValSans);
% hold on
% loglog(1./ValH, ValAvec);
% loglog(1./ValH, ValCea);
% legend('without postpro','WITH postpro','CEA degree P')


% valH = [0.5 0.25 0.125 0.0625 0.03125];
% A = [
%     0.69117 0.3256 0.68275 0.15757;
%     0.097293 0.046679 0.084009 0.011502;
%     0.00406 0.0032242 0.0013622 0.00037536;
%     0.00025464 0.00020996 3.3576e-05 1.1977e-05;
%     1.6187e-05 1.3437e-05 1.0293e-06 3.83e-07
%     ];
% 
% hold off
% loglog(1./valH, A(:,1)', 'b');
% hold on
% loglog(1./valH, A(:,3)', 'r');
% loglog(1./valH, A(:,2)', ':b');
% loglog(1./valH, A(:,4)', ':r');
% box on;
% grid on;
% xlabel('1/h')
% ylabel('error')
% legend('Without post-pro', 'With post-pro', 'Best approx. p=3', 'Best approx. p=4')
% 


% writeFieldDG(dofm, mesh, solP, "mySol.pos", "mySol");
% system('gmsh mySol.pos');

% figure(1);
% subplot(1,3,1);
% hold off
% postProVizu2D_DG(mesh,real(solP), 'Projected solution');
% subplot(1,3,2);
% hold off
% postProVizu2D_DG(mesh,real(solA), 'Numerical solution');
% subplot(1,3,3);
% hold off
% postProVizu2D_DG(mesh,real(solP-solA), 'Error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% fprintf('Solver  : gmres A\n');
% [solA,~,~,iterA]         = gmres(matA,rhsA,size(matA,1),resTol,size(matA,1));
% errorL2IterA             = computeNormError2D_DG(mesh, dofm, solA, solRef);
% fprintf('Solver  : bicgstab A\n');
% [solA,~,~,iterBiCGStabA] = bicgstab(matA,rhsA,resTol,size(matA,1));
% errorL2BiCGStabA         = computeNormError2D_DG(mesh, dofm, solA, solRef);
% fprintf('Solver  : conjgradn A\n');
% [solA,~,~,iterCGNA]      = conjgradn(matA,rhsA,resTol,size(matA,1));
% errorL2CGNA              = computeNormError2D_DG(mesh, dofm, solA, solRef);
% fprintf('Solver  : relaxation A\n');
% [solA,~,~,iterJacobiA]   = jacobi(matA,rhsA,resTol,size(matA,1),0.5);
% errorL2JacobiA           = computeNormError2D_DG(mesh, dofm, solA, solRef);

% matA11 = matA(1:3*dofm.numDofTRI, 1:3*dofm.numDofTRI);
% matA12 = matA(1:3*dofm.numDofTRI, 3*dofm.numDofTRI+1:end);
% matA21 = matA(3*dofm.numDofTRI+1:end, 1:3*dofm.numDofTRI);
% matA22 = matA(3*dofm.numDofTRI+1:end, 3*dofm.numDofTRI+1:end);
% rhsA1 = rhsA(1:3*dofm.numDofTRI);
% rhsA2 = rhsA(3*dofm.numDofTRI+1:end);
% matS = matA22 - matA21*(matA11\matA12);
% rhsS = rhsA2 - matA21*(matA11\rhsA1);

% fprintf('Solver  : gmres S\n');
% [solS,~,~,iterS]         = gmres(matS,rhsS,size(matS,1),resTol,size(matS,1));
% solA = matA11\(rhsA1 - matA12*solS);
% errorL2IterS             = computeNormError2D_DG(mesh, dofm, solA, solRef);
% fprintf('Solver  : bicgstab S\n');
% [solS,~,~,iterBiCGStabS] = bicgstab(matS,rhsS,resTol,size(matS,1));
% solA = matA11\(rhsA1 - matA12*solS);
% errorL2BiCGStabS         = computeNormError2D_DG(mesh, dofm, solA, solRef);
% fprintf('Solver  : conjgradn S\n');
% %[solS,flag,~,iterCGNS]      = conjgradn(matS,rhsS,resTol,10*size(matS,1));
% %solA = matA11\(rhsA1 - matA12*solS);
% %errorL2CGNS              = computeNormError2D_DG(mesh, dofm, solA, solRef);
% [solS,~,~,iterCGNS]      = conjgradnResPhy(matS,rhsS,resTol,10*size(matS,1),matA11,matA12,matA21,matA22,rhsA1,rhsA2);
% solA = matA11\(rhsA1 - matA12*solS);
% errorL2CGNS              = computeNormError2D_DG(mesh, dofm, solA, solRef);
% fprintf('Solver  : relaxation S\n');
% [solS,~,~,iterJacobiS]   = jacobi(matS,rhsS,resTol,size(matS,1),0.5);
% solA = matA11\(rhsA1 - matA12*solS);
% errorL2JacobiS           = computeNormError2D_DG(mesh, dofm, solA, solRef);

% fprintf('Solver  : richardson S\n');
% [solS,~,~,iterRichS]   = richardsonResPhy(matS,rhsS,resTol,10*size(matS,1),0.5,matA11,matA12,matA21,matA22,rhsA1,rhsA2);
% solA = matA11\(rhsA1 - matA12*solS);
% errorL2RichS           = computeNormError2D_DG(mesh, dofm, solA, solRef);

% fprintf('Solver  : eigenval A\n');
% [eigenvecA,eigenvalA] = eigs(matA,size(matA,1));
% eigenvalA = diag(eigenvalA);
% [eigenvecAA,eigenvalAA] = eigs(matA'*matA,size(matA,1));
% eigenvalAA = diag(eigenvalAA);

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
% disp(['    IterRich           ' num2str(iterRichS)]);
% disp(['    Final L2-Error     ' num2str(errorL2RichS)]);
% disp(['---------------------------------------------------------']);

% disp(['\text{HDG-1}(\tau=1) & full & ' ...
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

% disp(['\text{HDG-1}(\tau=i) & red & ' ...
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

