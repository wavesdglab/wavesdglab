%close all;
clear all;

global k

% hPower = [-1 -1.5 -2 -2.5 -3 -3.5];
% hList = 2.^hPower;
% for ITER = 1:size(hPower,2)
% h = hList(ITER);

% FREE SPACE
% k = 10*pi;
% h = 1/16;
% degree = 3;
% mesh = benchmark2D('open',h);

% WAVEGUIDE
% k = 6*pi;
% h = 1/8;
% degree = 3;
% mesh = benchmark2D('waveguide',h);

% CAVITY
k = 5.125*sqrt(2)*pi;
h = 1/16;
degree = 3;
mesh = benchmark2D('cavity',h);

% Define parameters
tol = 1e-4;
tau = 1;

% Build mesh and dofManager
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method HDG']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_HDG(mesh, dofm, tau);
[errorL2] = computeNormError2D_DG(mesh, dofm, solA);

[solP, sysP] = computeSolProjL2_2D(mesh, dofm);
[errorProjL2] = computeNormError2D_DG(mesh, dofm, solP);

[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solA);
[errorPostL2] = computeNormError2D_DG(mesh, dofmPost, solApost);

[solPpost, sysPpost] = computeSolProjL2_2D(mesh, dofmPost);
[errorProjPostL2] = computeNormError2D_DG(mesh, dofmPost, solPpost);

% disp([num2str(k) ' ' num2str(h) ' ' num2str(degree) ' ' num2str(Dlambda) ' ' num2str(errorL2) ' ' num2str(errorProjL2) ' ' num2str(errorPostL2) ' ' num2str(errorProjPostL2)]);

disp(['    L2-Error (numSol)   ' num2str(errorL2, '%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2, '%1.2e')]);
disp(['    L2-Error (numPost)  ' num2str(errorPostL2, '%1.2e')]);
disp(['    L2-Error (projPost) ' num2str(errorProjPostL2, '%1.2e')]);
disp('---------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(['--- CALL eigsA']);
% [eigenvecA,eigenvalA] = eigs(sysA.matA,size(sysA.matA,1));
% eigenvalA = diag(eigenvalA);
% 
% disp(['A : Size                ' num2str(size(sysA.matA,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecA))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecA))]);
% disp(['    Cond(A)             ' num2str(condest(sysA.matA))]);
% disp(['---------------------------------------------------------']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(['--- CALL eigsS']);
% [eigenvecS,eigenvalS]      = eigs(sysA.matS,size(sysA.matS,1));
% eigenvalS                  = diag(eigenvalS);
% 
% disp(['S : Size                ' num2str(size(sysA.matS,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecS))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecS))]);
% disp(['    Cond(S)             ' num2str(condest(sysA.matS), '%1.2e')]);
% disp(['    Cond(SS)            ' num2str(condest(sysA.matS'*sysA.matS), '%1.2e')]);
% disp(['---------------------------------------------------------']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

color1 = [0, 0.4470, 0.7410];      % blue
color2 = [0.8500, 0.3250, 0.0980]; % red
color3 = [0.9290, 0.6940, 0.1250]; % yellow
color4 = [0.4940, 0.1840, 0.5560]; % magenta
color5 = [0.4660, 0.6740, 0.1880]; % green
color6 = [0.3010, 0.7450, 0.9330]; % cyan
color7 = [0.6350, 0.0780, 0.1840]; % brown

disp(['--- CALL conjgradn']);
maxit = 1000;
itout = 10;
[resVec, resRedVec, resPhyVec, error, errorPP, iter, flag] = solverCGNreduDG(mesh, dofm, sysA, 1e-10, maxit, itout);

figure(3);
hold off
semilogy(0:itout:maxit,resVec,'Color',color4,'DisplayName','Relative residual CGN');
hold on
semilogy(0:itout:maxit,resPhyVec,'Color',color2,'DisplayName','Relative residual (Phy)');
semilogy(0:itout:maxit,resRedVec,'Color',color1,'DisplayName','Relative residual (Red)');
semilogy(0:itout:maxit,error,'Color',color3,'DisplayName','Relative L2-error');
semilogy(0:itout:maxit,errorPP,'Color',color5,'DisplayName','Relative L2-error with PostPro');
plot([0 maxit],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (Ref)');
plot([0 maxit],[errorPostL2 errorPostL2],'k--','DisplayName','Relative L2-error with PostPro (Ref)');
box on;
grid on;
legend('Location','southwest');
xlim([0 maxit]);
ylim([1e-5 1e0]);
xlabel('Iteration');
ylabel('Value');
set(gcf, 'PaperUnits', 'points','PaperPosition', [0 0 400 300]);
print(['~/Desktop/HistoConvCGN-BenchCavity-HDG.eps'],'-depsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% writeFieldDG(dofm, mesh, solP, "mySol.pos", "mySol");
% system('gmsh mySol.pos');



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

% disp(['\text{HDG}(\tau=1) & full & ' ...
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

% disp(['\text{HDG}(\tau=i) & red & ' ...
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

% end