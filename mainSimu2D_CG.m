%close all;
clear all;

headers2D();
global k

% hPower = [-1 -1.5 -2 -2.5 -3 -3.5 -4 -4.5];
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

% Build mesh and dofManager
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method CG']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_CG(mesh, dofm);
[errorL2, errorH1] = computeNormError2D_CG(mesh, dofm, solA);

[solP, sysP] = computeSolProjL2_2D_CG(mesh, dofm);
[errorProjL2, errorProjH1] = computeNormError2D_CG(mesh, dofm, solP);

% disp([num2str(k) ' ' num2str(h) ' ' num2str(degree) ' ' num2str(Dlambda) ' ' num2str(errorL2) ' ' num2str(errorProjL2)]);

disp(['    L2-Error (numSol)   ' num2str(errorL2, '%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2, '%1.2e')]);
disp(['---------------------------------------------------------']);

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

% color1 = [0, 0.4470, 0.7410];      % blue
% color2 = [0.8500, 0.3250, 0.0980]; % red
% color3 = [0.9290, 0.6940, 0.1250]; % yellow
% color4 = [0.4940, 0.1840, 0.5560]; % magenta
% color5 = [0.4660, 0.6740, 0.1880]; % green
% color6 = [0.3010, 0.7450, 0.9330]; % cyan
% color7 = [0.6350, 0.0780, 0.1840]; % brown
% 
% disp(['--- CALL conjgradn']);
% maxit = 1000;
% itout = 10;
% [resVec, resRedVec, resPhyVec, error, iter, flag] = solverCGNreduCG(mesh, dofm, sysA, 1e-10, maxit, itout);
% 
% figure(1);
% hold off
% semilogy(0:itout:maxit,resVec,'Color',color4,'DisplayName','Relative residual CGN');
% hold on
% semilogy(0:itout:maxit,resPhyVec,'Color',color2,'DisplayName','Relative residual (Phy)');
% semilogy(0:itout:maxit,resRedVec,'Color',color1,'DisplayName','Relative residual (Red)');
% semilogy(0:itout:maxit,error,'Color',color3,'DisplayName','Relative L2-error');
% plot([0 maxit],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (Ref)');
% box on;
% grid on;
% legend('Location','northeast');
% xlim([0 maxit]);
% ylim([1e-5 1e0]);
% xlabel('Iteration');
% ylabel('Value');
% set(gcf, 'PaperUnits', 'points','PaperPosition', [0 0 400 300]);
% print(['~/Desktop/HistoConvCGN-BenchCavity-CG.eps'],'-depsc');
% legend('Location','northeast');

% disp(['--- CALL gmres']);
% [solAiter,~,~,iterA]         = gmres(sysA.matA,sysA.rhsA,size(sysA.matA,1),resTol,size(sysA.matA,1));
% errorL2IterA                 = computeNormError2D_CG(mesh, dofm, solAiter, solA);
% disp(['--- CALL bicgstab']);
% [solAiter,~,~,iterBiCGStabA] = bicgstab(sysA.matA,sysA.rhsA,resTol,100*size(sysA.matA,1));
% errorL2BiCGStabA             = computeNormError2D_CG(mesh, dofm, solAiter, solA);
% disp(['--- CALL conjgradn']);
% [solAiter,~,~,iterCGNA]      = solverCGN(sysA,resTol);
% errorL2CGNA                  = computeNormError2D_CG(mesh, dofm, solAiter, solA);
% disp(['--- CALL jacobi']);
% [solAiter,~,~,iterJacobiA]   = jacobi(sysA.matA,sysA.rhsA,resTol,size(sysA.matA,1),0.5);
% errorL2JacobiA               = computeNormError2D_CG(mesh, dofm, solAiter, solA);
% disp(['--- CALL eigenvalues']);
% 
% disp(['---------------------------------------------------------']);
% disp(['    IterGmres           ' num2str(iterA(2))]);
% disp(['    Final L2-Error      ' num2str(errorL2IterA)]);
% disp(['    IterBiCGS           ' num2str(iterBiCGStabA)]);
% disp(['    Final L2-Error      ' num2str(errorL2BiCGStabA)]);
% disp(['    IterCGN             ' num2str(iterCGNA)]);
% disp(['    Final L2-Error      ' num2str(errorL2CGNA)]);
% disp(['    IterRelax           ' num2str(iterJacobiA)]);
% disp(['    Final L2-Error      ' num2str(errorL2JacobiA)]);
% disp(['---------------------------------------------------------']);
% 
% disp(['\text{CG} & & ' ...
%     num2str(errorL2,'%.1e') ' & ' ...
%     num2str(size(sysA.matA,1)) ' & ' ...
%     num2str(rank(eigenvecA)) ' & ' ...
%     num2str(cond(eigenvecA),'%.1e') ' & ' ...
%     num2str(condest(sysA.matA),'%.1e') ' & ' ...
%     num2str(iterA(2)) ' & ' ...
%     num2str(errorL2IterA,'%.1e') ' & ' ...
%     num2str(iterBiCGStabA) ' & ' ...
%     num2str(errorL2BiCGStabA,'%.1e') ' & ' ...
%     num2str(iterCGNA) ' & ' ...
%     num2str(errorL2CGNA,'%.1e') ' & ' ...
%     num2str(iterJacobiA) ' & ' ...
%     num2str(errorL2JacobiA,'%.1e') ' \\'
%     ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

writeField_CG(dofm, mesh, solA, "solution.pos", "solution");
system('gmsh solution.pos');

% end
