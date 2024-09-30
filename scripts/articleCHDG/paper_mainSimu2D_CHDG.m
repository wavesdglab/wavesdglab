clear all;
%close all;

global omega eta k c rho eta1 eta2 k1 k2 c1 c2 rho1 rho2 h
   
% Setup benchmark and parameters
benchmark = 'open_heterogeneous';
switch benchmark
    case 'open_heterogeneous'
        omega = 15*pi; %15*pi;
        h = 1/16;  % to be set directly in open.geo 
        tol = 1e-10; maxit = 1000; itout = 50;
        rho1 = 1;
        c1 = 1;  
        rho2 = 2;
        c2 = 1/2; 
        eta1 = rho1 * c1;
        eta2 = rho2 * c2;
        k1 = omega / c1;
        k2 = omega / c2;
    case 'disk_heterogeneous'
        omega = 10*pi; %10*pi;
        h = 0.055; % to be set directly in disk.geo
        tol = 1e-10; maxit = 1000; itout = 100;
        rho1 = 1;
        c1 = 1; 
        rho2 = 1;
        c2 = 1;
        eta1 = rho1 * c1;
        eta2 = rho2 * c2;
        k1 = omega / c1;
        k2 = omega / c2;
end
degree = 3;
BASIS = 0;
PREC = 1;
% A=1 and B=1 for 0th order, A=2 and B=2 for 2nd order
A = 2;              % order of numerical fluxes
B = 2;              % order of transmission variables

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);
setParameters(mesh, benchmark);
Dlambda = 2*pi/k(1) * (sqrt(dofm.numDofTRI) - 1);

%%
% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CHDG']);
disp(['---------------------------------------------------------']);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

% Compute numerical solution
% [solA, sysA] = computeSolNum2D_CHDG_ALL(mesh, dofm, PREC, A, B);
[solA, sysA] = computeSolNum2D_CHDG_heterogeneous_upw(mesh, dofm, PREC);
% [solA, sysA] = computeSolNum2D_CHDG_heterogeneous_mean(mesh, dofm, 0, 1);
% [solA, sysA] = computeSolNum2D_HDG_ALL(mesh, dofm, BASIS, PREC);   % TO BE COMPARED

%%
% Compute numerical error
errorL2_A = computeNormError2D_DG_ALL(mesh, dofm, solA);
% errorL2_B = computeNormError2D_DG_ALL(mesh, dofm, solB)
% errorL2_C = computeNormError2D_DG_ALL(mesh, dofm, solC)
errorL2 = errorL2_A

% Compute projection solution
solP = computeSolProjL2_2D_DG(mesh, dofm);
% errorProjL2 = computeNormError2D_DG_ALL(mesh, dofm, solP);
% 
% disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
% disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
% disp('---------------------------------------------------------');

%%
% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------
%
writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
writeField2D(dofm, mesh, solA(1:mesh.numTri*3*dofm.numDofPerTRI)-solP, 'output/errNum.pos', "errNum");
system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

%%
% -------------------------------------------------------------------------
% Compute eigenvalues/eigenvectors
% -------------------------------------------------------------------------

% mat = sysA.matPinv*sysA.matS;
% [~, eigenval] = eigs(mat,size(mat,1));
% 
% eigenval = 1 - diag(eigenval);
% 
% 1 - max(abs(eigenval))
% 
% fprintf('Spectral radius = %.16f\n', max(abs(eigenval)));

%%
% -------------------------------------------------------------------------
% Compute iterative solution
% -------------------------------------------------------------------------

% solver = 'GMRES';
% switch solver
%     case 'Rich'
%         alpha = 1;
%         [resRedVec, resPhyVec, errorVec] = solverRichardsonRedu_DG(mesh, dofm, sysA, tol, maxit, itout, alpha, @computeNormError2D_DG_ALL);
%     case 'CGNR'
%         [resRedVec, resPhyVec, errorVec] = solverCGNRredu_DG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_DG_ALL);
%     case 'GMRES'
%         [resRedVec, resPhyVec, errorVec] = solverGMRESredu_DG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_DG_ALL);
% 
% end
 
% alpha = 1;
% [CHDG_resRedVec_Rich, CHDG_resPhyVec_Rich, CHDG_errorVec_Rich] = solverRichardsonRedu_DG(mesh, dofm, sysA, tol, maxit, itout, alpha, @computeNormError2D_DG_ALL);
% [CHDG_resRedVec_CGNR, CHDG_resPhyVec_CGNR, CHDG_errorVec_CGNR] = solverCGNRredu_DG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_DG_ALL);
% [CHDG_resRedVec_GMRES, CHDG_resPhyVec_GMRES, CHDG_errorVec_GMRES] = solverGMRESredu_DG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_DG_ALL);
% [HDG_resRedVec_CGNR, HDG_resPhyVec_CGNR, HDG_errorVec_CGNR] = solverCGNRredu_DG(mesh, dofm, sysB, tol, maxit, itout, @computeNormError2D_DG_ALL);
% [HDG_resRedVec_GMRES, HDG_resPhyVec_GMRES, HDG_errorVec_GMRES] = solverGMRESredu_DG(mesh, dofm, sysB, tol, maxit, itout, @computeNormError2D_DG_ALL);


% figure();
% hold off
% % semilogy(0:itout:maxit,resPhyVec,'r','DisplayName','Relative residual (Phy)');
% % semilogy(0:itout:maxit,resRedVec,'b','DisplayName','Relative residual (Red)');
% % semilogy(0:itout:maxit,errorVec(1,:),'r-o','DisplayName','Relative L2-error (iterative - 1st order)');
% % hold on
% % semilogy(0:itout:maxit,errorVec(2,:),'g-o','DisplayName','Relative L2-error (iterative - 2nd order)');
% % hold on
% semilogy(0:itout:maxit,errorVec(3,:),'r-x','DisplayName','Relative L2-error (iterative - 1st order) - precond.');
% hold on
% semilogy(0:itout:maxit,errorVec(4,:),'g-x','DisplayName','Relative L2-error (iterative - 2nd order) - precond.');
% hold on
% plot([0 maxit],[errorL2_1 errorL2_1],'k--','DisplayName','Relative L2-error (direct - 1st order)');
% hold on
% plot([0 maxit],[errorL2_2 errorL2_2],'k-.','DisplayName','Relative L2-error (direct - 2nd order)');
% % plot([0 maxit],[errorProjL2 errorProjL2],'k:','DisplayName','Relative L2-error (projection)');
% box on;
% grid on;
% % legend('Location','southwest');
% legend('Iterative - 1st order','Iterative - 2nd order', 'Iterative - 1st order - precond.', 'Iterative - 2nd order - precond.', 'Direct - 1st order', 'Direct - 2nd order');
% title(['CHDG - ' benchmark ' - ' solver ' - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)])
% xlim([0 maxit]);
% ylim([0.005 1]);
% xlabel('Iteration');
% ylabel('Relative L^2-error');

% figure();
% hold off
% semilogy(0:itout:maxit,errorVec_1,'r-x','DisplayName','Relative L2-error (iterative - 1st order) - precond.');
% hold on
% semilogy(0:itout:maxit,errorVec_2,'g-x','DisplayName','Relative L2-error (iterative - 2nd order) - precond.');
% hold on
% plot([0 maxit],[errorL2_1 errorL2_1],'k--','DisplayName','Relative L2-error (direct - 1st order)');
% hold on
% plot([0 maxit],[errorL2_2 errorL2_2],'k-.','DisplayName','Relative L2-error (direct - 2nd order)');
% box on;
% grid on;
% legend('Location','southwest');
% legend('Iterative - 1st order - precond.', 'Iterative - 2nd order - precond.', 'Direct - 1st order', 'Direct - 2nd order');
% title(['CHDG - ' benchmark ' - ' solver ' - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)])
% xlim([0 maxit]);
% ylim([0.005 1]);
% xlabel('Iteration');
% ylabel('Relative L^2-error');
 
% figure();
% hold off
% semilogy(0:itout:maxit,HDG_errorVec_CGNR,'-or','MarkerFaceColor','w','DisplayName','HDG - CGN');
% hold on
% semilogy(0:itout:maxit,HDG_errorVec_GMRES,'-or','MarkerFaceColor','r','DisplayName','HDG - GMRES');
% semilogy(0:itout:maxit,CHDG_errorVec_CGNR,'-ob','MarkerFaceColor','w','DisplayName','CHDG - CGN');
% hold on
% semilogy(0:itout:maxit,CHDG_errorVec_GMRES,'-ob','MarkerFaceColor','b','DisplayName','CHDG - GMRES');
% semilogy(0:itout:maxit,CHDG_errorVec_Rich,'-xb','DisplayName','CHDG - Fixed-point');
% semilogy([0 maxit],[errorL2 errorL2],'k--','DisplayName','Direct solver');
% box on;
% grid on;
% legend('Location','southwest');
% legend('Location','northoutside','NumColumns',4);
% %title(['CHDG - ' benchmark ' - \omega=' num2str(omega) ' - h=' num2str(h) ' - P=' num2str(degree)])
% xlim([0 maxit]);
% ylim([0.05 1]);
% xlabel('Iteration');
% ylabel('Relative error');