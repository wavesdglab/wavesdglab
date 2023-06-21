clear all;
%close all;

global k

% Setup benchmark and parameters
benchmark = 'cavity';
switch benchmark
    case 'open'
        k = 30*pi; %15*pi;
        h = 1/34; %1/16;
        tol = 1e-10; maxit = 1000; itout = 50;
    case 'cavity'
        k = 7.01*sqrt(2)*pi; %7.1*sqrt(2)*pi;
        h = 1/15; %1/10;
        tol = 1e-10; maxit = 2000; itout = 100;
    case 'waveguide'
        k = 12*pi; %6*pi
        h = 1/17; %1/8
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 3;
tau = 1;
BASIS = 1;
PREC = 0;

% Build mesh and DOF manager
mesh = benchmark2D(benchmark,h);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CHDG']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['    tau                 ' num2str(tau)]);
disp(['---------------------------------------------------------']);

% Compute numerical solution/error
[solA, sysA] = computeSolNum2D_CHDG(mesh, dofm, tau, BASIS, PREC);
errorL2 = computeNormError2D_DG(mesh, dofm, solA);

% Compute projection solution/error
solP = computeSolProjL2_2D_DG(mesh, dofm);
errorProjL2 = computeNormError2D_DG(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
disp('---------------------------------------------------------');

% disp([num2str(k) ' ' num2str(h) ' ' num2str(degree) ' ' num2str(Dlambda) ' ' num2str(errorL2) ' ' num2str(errorProjL2)]);

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
writeField2D(dofm, mesh, solA-solP, 'output/errNum.pos', "errNum");
system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

% -------------------------------------------------------------------------
% Compute eigenvalues/eigenvectors
% -------------------------------------------------------------------------

% mat = sysA.matPinv*sysA.matS;
% [eigenvec, eigenval] = eigs(mat,size(mat,1));
% eigenval = diag(eigenval);
% rankEigenVec = rank(eigenvec);
% condEigenVec = cond(eigenvec);
% 
% disp(['    Size                ' num2str(size(mat,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rankEigenVec)]);
% disp(['    Cond(eigenvectors)  ' num2str(condEigenVec)]);
% 
% figure;
% hold off; scatter(real(eigenval),imag(eigenval));
% %hold on; plot(fovals(mat,100));
% hold on; plot(cos(0:0.01:2*pi)+1,sin(0:0.01:2*pi),'k');
% grid on; box on;
% set(gcf, 'PaperUnits', 'points','PaperPosition', [0 0 500 500]);
% print(['output/Eigenvalues-' benchmark '-CHDG.png'],'-dpng');
% 
% % Compute condition number
% condestMat = condest(mat);
% 
% disp(['    Cond(mat)           ' num2str(condestMat,'%1.2e')]);
% disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Compute iterative solution
% -------------------------------------------------------------------------

% solver = 'CGNR';
% switch solver
%     case 'Rich'
%         alpha = 1;
%         [resRedVec, resPhyVec, errorVec] = solverRichardson_DG(mesh, dofm, sysA, tol, maxit, itout, alpha, @computeNormError2D_DG);
%     case 'CGNR'
%         [resRedVec, resPhyVec, errorVec] = solverCGNRredu_DG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_DG);
%     case 'GMRES'
%         [resRedVec, resPhyVec, errorVec] = solverGMRESredu_DG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_DG);
% end
% 
% figure;
% hold off
% semilogy(0:itout:maxit,resPhyVec,'r','DisplayName','Relative residual (Phy)');
% hold on
% semilogy(0:itout:maxit,resRedVec,'b','DisplayName','Relative residual (Red)');
% semilogy(0:itout:maxit,errorVec,'k','DisplayName','Relative L2-error (iterative)');
% plot([0 maxit],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (direct)');
% plot([0 maxit],[errorProjL2 errorProjL2],'k:','DisplayName','Relative L2-error (projection)');
% box on;
% grid on;
% legend('Location','southwest');
% title(['CHDG - ' benchmark ' - ' solver ' - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)])
% xlim([0 maxit]);
% ylim([0.005 1]);
% xlabel('Iteration');
% ylabel('Value');