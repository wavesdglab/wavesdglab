%close all;
clear all;

headers2D;
global k

% hPower = [-1 -1.5 -2 -2.5 -3 -3.5];
% hList = 2.^hPower;
% for ITER = 1:size(hPower,2)
% h = hList(ITER);

% EIGEN OPEN
% k = 5*pi;
% h = 1/8;
% degree = 3;
% mesh = benchmark2D('open',h);

% EIGEN WAVEGUIDE
% k = 5*pi;
% h = 1/8;
% degree = 3;
% mesh = benchmark2D('waveguide',h);

% EIGEN CLOSE
% k = 5*pi;
% h = 1/8;
% degree = 3;
% mesh = benchmark2D('cavity',h);

% BENCH FREE SPACE
k = 10*pi;
h = 1/8;
degree = 3;
mesh = benchmark2D('open',h);

% BENCH WAVEGUIDE
% k = 6*pi;
% h = 1/8;
% degree = 3;
% mesh = benchmark2D('waveguide',h);

% BENCH CAVITY
% k = 5.125*sqrt(2)*pi;
% h = 1/20;
% degree = 3;
% mesh = benchmark2D('cavity',h);

% Define parameters
tol = 1e-4;
tau = 0.1;

% Build mesh and dofManager
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method UDG-1']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['    tau                 ' num2str(tau)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_UDG1(mesh, dofm, tau);
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

% writeFieldDG(dofm, mesh, solP, "mySol.pos", "mySol");
% system('gmsh mySol.pos');

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

disp(['--- CALL eigsS']);
[eigenvecS,eigenvalS]      = eigs(sysA.matS,size(sysA.matS,1));
eigenvalS                  = diag(eigenvalS);

% disp(['S : Size                ' num2str(size(sysA.matS,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecS))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecS))]);
% disp(['    Cond(S)             ' num2str(condest(sysA.matS), '%1.2e')]);
% disp(['    Cond(SS)            ' num2str(condest(sysA.matS'*sysA.matS), '%1.2e')]);
% disp(['---------------------------------------------------------']);
% 
figure(3);
hold off
scatter(real(eigenvalS),imag(eigenvalS));
%hold on
%plot(fovals(sysA.matS,100),'-b','DisplayName','Numerical range');
grid on; box on;
%title(['Eigenvalues S  —  Min real eigenvalues: ' num2str(min(real(eigenvalS)))]);
%legend();
axis([-0.02 0.14 -0.08 0.08]);
set(gcf, 'PaperUnits', 'points','PaperPosition', [0 0 500 500]);
print(['~/Desktop/Eigenvalues-BenchCavity-UDG1.eps'],'-depsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% color1 = [0, 0.4470, 0.7410];      % blue
% color2 = [0.8500, 0.3250, 0.0980]; % red
% color3 = [0.9290, 0.6940, 0.1250]; % yellow
% color4 = [0.4940, 0.1840, 0.5560]; % magenta
% color5 = [0.4660, 0.6740, 0.1880]; % green
% color6 = [0.3010, 0.7450, 0.9330]; % cyan
% color7 = [0.6350, 0.0780, 0.1840]; % brown
% color8 = [0.3, 0.3, 0.3]; % gray
% color9 = [0.0, 0.0, 0.0]; % gray
% 
% disp(['--- CALL solverCGNreduDG']);
% maxit = 1000;
% itout = 20;
% [resVec, resRedVec, resPhyVec, error, errorPP, iter, flag] = solverCGNreduDG(mesh, dofm, sysA, 1e-10, maxit, itout);
% 
% figure(4);
% hold off
% semilogy(0:itout:maxit,resPhyVec,'o-','Color',color2,'DisplayName','Relative residual - Physical system');
% hold on
% semilogy(0:itout:maxit,resRedVec,'x-','Color',color1,'DisplayName','Relative residual - Reduced system');
% semilogy(0:itout:maxit,error,'o-','Color',color3,'DisplayName','Relative numerical error');
% semilogy(0:itout:maxit,errorPP,'x-','Color',color5,'DisplayName','Relative numerical error with post-processing');
% plot([0 maxit],[errorL2 errorL2],'--','Color',color8,'DisplayName','Relative projection error');
% plot([0 maxit],[errorPostL2 errorPostL2],'-','Color',color9,'DisplayName','Relative projection error with post-processing');
% box on;
% grid on;
% legend('Location','southwest');
% xlim([0 maxit]);
% ylim([1e-5 1e0]);
% xlabel('Iteration');
% ylabel('Value');
% set(gcf, 'PaperUnits', 'points','PaperPosition', [0 0 400 300]);
% print(['~/Desktop/HistoConvCGN-BenchWaveguide-UDGi.eps'],'-depsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(['--- CALL solverRichardsonDG']);
% maxit = 100000;
% itout =   5000;
% 
% alpha1 = 1;
% alpha2 = 0.8;
% [resRedVec1, resPhyVec1, error1, errorPP1] = solverRichardsonDG(mesh, dofm, sysA, 1e-10, maxit, itout, alpha1);
% [resRedVec2, resPhyVec2, error2, errorPP2] = solverRichardsonDG(mesh, dofm, sysA, 1e-10, maxit, itout, alpha2);
% 
% figure(3);
% hold off
% semilogy(0:itout:maxit,resPhyVec1,'-x','Color',color2,'DisplayName','Relative residual (Phy) — \alpha = 1');
% hold on
% semilogy(0:itout:maxit,resPhyVec2,'-o','Color',color2,'DisplayName','Relative residual (Phy) — \alpha = 0.8');
% semilogy(0:itout:maxit,resRedVec1,'-x','Color',color1,'DisplayName','Relative residual (Red) — \alpha = 1');
% semilogy(0:itout:maxit,resRedVec2,'-o','Color',color1,'DisplayName','Relative residual (Red) — \alpha = 0.8');
% semilogy(0:itout:maxit,error1,'-x','Color',color3,'DisplayName','Relative L2-error — \alpha = 1');
% semilogy(0:itout:maxit,error2,'-o','Color',color3,'DisplayName','Relative L2-error — \alpha = 0.8');
% semilogy(0:itout:maxit,errorPP1,'-x','Color',color5,'DisplayName','Relative L2-error with PostPro — \alpha = 1');
% semilogy(0:itout:maxit,errorPP2,'-o','Color',color5,'DisplayName','Relative L2-error with PostPro — \alpha = 0.8');
% plot([0 maxit],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (Ref)');
% plot([0 maxit],[errorPostL2 errorPostL2],'k--','DisplayName','Relative L2-error with PostPro (Ref)');
% box on;
% grid on;
% legend('Location','southwest');
% xlim([0 maxit]);
% ylim([1e-5 1e0]);
% xlabel('Iteration');
% ylabel('Value');
% set(gcf, 'PaperUnits', 'points','PaperPosition', [0 0 400 300]);
% print(['~/Desktop/HistoConvRichardson-BenchWaveguide-UDG.eps'],'-depsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('Solver  : gmres A\n');
% [solGiter,~,~,iterS] = gmres(sysA.matS,sysA.rhsS,size(sysA.matS,1),resTol,size(sysA.matS,1));
% solIiter = sysA.matIIinv*(sysA.rhsI-sysA.matIG*solGiter);
% solAiter = [ solIiter ; solGiter ];
% errorL2IterA = computeNormError2D_DG(mesh, dofm, solAiter, solA);

% fprintf('Solver  : bicgstab A\n');
% [solSiter,~,~,iterBiCGStabA] = bicgstab(matA,rhsA,resTol,size(matA,1));
% errorL2BiCGStabA         = computeNormError2D_DG(mesh, dofm, solAiter, solA);
% fprintf('Solver  : conjgradn A\n');
% [solSiter,~,~,iterCGNA]      = conjgradn(matA,rhsA,resTol,size(matA,1));
% errorL2CGNA              = computeNormError2D_DG(mesh, dofm, solAiter, solA);
% fprintf('Solver  : relaxation A\n');
% [solSiter,~,~,iterJacobiA]   = jacobi(matA,rhsA,resTol,size(matA,1),0.5);
% errorL2JacobiA           = computeNormError2D_DG(mesh, dofm, solAiter, solA);

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
% [eigenvecA,eigenvalA] = eigs(sysA.matA,size(sysA.matA,1));
% eigenvalA = diag(eigenvalA);
% 
% fprintf('Solver  : eigenval S\n');
% [eigenvecS,eigenvalS] = eigs(sysA.matS,size(sysA.matS,1));
% eigenvalS = diag(eigenvalS);
% 
% figure(1);
% hold off
% scatter(real(eigenvalA),imag(eigenvalA),'b','DisplayName','Eigenvalues');
% hold on
% %plot(fovals(sysA.matA,100),'-b','DisplayName','Numerical range');
% grid on; box on;
% title(['Eigenvalues A']);
% legend();
% axis([-0.05 0.25 -0.25 0.1]);
% 
% figure(2);
% hold off
% scatter(real(eigenvalS),imag(eigenvalS),'b','DisplayName','Eigenvalues');
% hold on
% %plot(fovals(sysA.matS,100),'-b','DisplayName','Numerical range');
% grid on; box on;
% title(['Eigenvalues S']);
% legend();
% axis([-0.05 0.25 -0.25 0.1]);

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

% end