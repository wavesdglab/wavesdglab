clear all;
%close all;

global omega eta k eta1 eta2 k1 k2 c1 c2 rho1 rho2 h

% Setup benchmark and parameters
benchmark = 'disk_heterogeneous';
switch benchmark
    case 'open'
        k = 15*pi; %15*pi;
        h = 1/16;  % 1/16
        tol = 1e-10; maxit = 200; itout = 50;  
    case 'cavity'
        k = 7.1*sqrt(2)*pi; %7.01*sqrt(2)*pi;
        h = 0.1; %1/10;
        tol = 1e-10; maxit = 2000; itout = 100;
    case 'waveguide'
        k = 6*pi; %6*pi
        h = 1/6; %1/8
        tol = 1e-10; maxit = 4000; itout = 200;
    case 'open_heterogeneous'
        omega = 15*pi; %15*pi;
        h = 1/16;  % 1/16
        tol = 1e-10; maxit = 1000; itout = 50;
        rho1 = 1;
        c1 = 1;  
        rho2 = 1;
        c2 = 1;  
        eta1 = rho1 * c1;
        eta2 = rho2 * c2;
        k1 = omega / c1;
        k2 = omega / c2;
    case 'cavity_heterogeneous'
        omega = 15*pi; 
        h = 1/16;
        tol = 1e-10; maxit = 1000; itout = 50;
        rho1 = 1;
        c1 = 1.3;
        rho2 = 1;
        c2 = 0.8;
        eta1 = rho1 * c1;
        eta2 = rho2 * c2;
        k1 = omega / c1;
        k2 = omega / c2;
    case 'waveguide_heterogeneous'
        omega = 6*pi; %6*pi
        h = 1/4; %1/8
        tol = 1e-10; maxit = 4000; itout = 200;
        rho1 = 1;
        c1 = 2;
        rho2 = 1;
        c2 = 1.2;
        eta1 = rho1 * c1;
        eta2 = rho2 * c2;
        k1 = omega / c1;
        k2 = omega / c2;
    case 'disk_heterogeneous'
        omega = 10*pi; 
        h = 0.0125;
        tol = 1e-10; maxit = 1000; itout = 100;
        rho1 = 1;
        c1 = 1; 
        rho2 = 1;
        c2 = 1/2;
        eta1 = rho1 * c1;
        eta2 = rho2 * c2;
        k1 = omega / c1;
        k2 = omega / c2;
end
degree = 3;
BASIS = 0;
PREC = 1;
A = 1;              % order of numerical fluxes
B = 1;              % order of transmission variables

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);
setParameters(mesh, benchmark);
Dlambda = 2*pi/k(1) * (sqrt(dofm.numDofTRI) - 1);

disp('Num');
mesh.numTri

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

% Compute numerical solution (Upwind)
[solA, sysA] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, BASIS, PREC);
% [solB, sysB] = computeSolNum2D_HDG_heterogeneous(mesh, dofm, BASIS, PREC);
% PREC = 0;
% [solC, sysC] = computeSolNum2D_DG_heterogeneous(mesh, dofm, PREC);

% Compute numerical solution (High-order)
% [solA, sysA] = computeSolNum2D_CHDG_ALL(mesh, dofm, PREC, A, B);
% [solB, sysB] = computeSolNum2D_HDG_ALL(mesh, dofm, BASIS, PREC);
% PREC = 0;
% [solC, sysC] = computeSolNum2D_DG_ALL(mesh, dofm, PREC);

%%
% Compute numerical error
errorL2_A = computeNormError2D_DG_ALL(mesh, dofm, solA);
% % errorL2_B = computeNormError2D_DG_ALL(mesh, dofm, solB)
% % errorL2_C = computeNormError2D_DG_ALL(mesh, dofm, solC)
errorL2 = errorL2_A;
% 
% Compute projection solution
solP = computeSolProjL2_2D_DG(mesh, dofm);
errorProjL2 = computeNormError2D_DG_ALL(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
disp('---------------------------------------------------------');

% figure()
% loglog(H,err(:,1),'-r','LineWidth',2);
% hold on
% loglog(H,err(:,2),'-b','LineWidth',2);
% loglog(H,err(:,3),'-g','LineWidth',2);
% title('L^2 error');
% legend('1^{st} order', '2^{nd} order', 'Projection', 'Location', 'northwest');
% grid on
% xlabel('h');
% ylabel('L^2 error');

% disp([num2str(k) ' ' num2str(h) ' ' num2str(degree) ' ' num2str(Dlambda) ' ' num2str(errorL2) ' ' num2str(errorProjL2)]);

%%
% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

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
% % eigenval = diag(eigenval);
% 
% eigenval = 1 - diag(eigenval);
% 
% 1 - max(abs(eigenval))

% A=0;
% B=0;
% m=0;
% n=0;
% for i=1:size(eigenval)
%     norm = abs(eigenval(i,1));
%     if norm>1
%         m=m+1;
%         A(m)=eigenval(i,1);
%     else
%         n=n+1;
%         B(n)=eigenval(i,1);
%     end
% end

% fprintf('Spectral radius = %.16f\n', max(abs(eigenval)));
% fprintf('Number of eigenvalues outside the unit circle = %i\n', m);
% fprintf('Number of eigenvalues = %i\n', i);

% rankEigenVec = rank(eigenvec);
% condEigenVec = cond(eigenvec);
% % 
% disp(['    Size                ' num2str(size(mat,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rankEigenVec)]);
% disp(['    Cond(eigenvectors)  ' num2str(condEigenVec)]);

% figure;
% hold off; scatter(real(eigenval),imag(eigenval),"blue");
% axis equal
% % hold on; plot(fovals(mat,100));
% hold on; %plot(cos(0:0.01:2*pi)+1,sin(0:0.01:2*pi),'k');
% plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k');
% grid on; box on;
% set(gcf, 'PaperUnits', 'points','PaperPosition', [0 0 500 500]);
% print(['output/Eigenvalues-' benchmark '-CHDG.png'],'-dpng');

% % Plot spectrum
% if A==0
%     figure;
%     hold off;scatter(real(eigenval),imag(eigenval));
%     % hold on; plot(fovals(mat,100));
%     hold on; plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k');
%     grid on; box on; axis equal;
% else
%     % Plot spectrum
%     figure;
%     hold off;
%     plot(real(A),imag(A),'xr');
%     hold on
%     plot(real(B),imag(B),'ob');
%     plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k');
%     grid on; box on; axis equal;
%     legend('|\lambda|>1','|\lambda|<1','FontSize',12);
% end

% % Compute condition number
% condestMat = condest(mat);
% 
% disp(['    Cond(mat)           ' num2str(condestMat,'%1.2e')]);
% disp(['---------------------------------------------------------']);

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
% [DG_resVec_CGNR, DG_errorVec_CGNR] = solverCGNR(mesh, dofm, sysC, tol, maxit, itout, @computeNormError2D_DG_ALL);
% [DG_resVec_GMRES, DG_errorVec_GMRES] = solverGMRES(mesh, dofm, sysC, tol, maxit, itout, @computeNormError2D_DG_ALL);

% figure(26);
% plot(rr,errorVec,'r-o');
% box on;
% grid on;
% % legend('Location','southwest');
% % legend('Iterative - 1st order - precond.', 'Iterative - 2nd order - precond.', 'Direct - 1st order', 'Direct - 2nd order');
% title(['CHDG - ' benchmark ' - ' solver ' - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)])
% xlim([gamma_min gamma_max]);
% ylim([0.005 0.1]);
% xlabel('\gamma');
% ylabel('Relative L^2-error');

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
 
% figure;
% hold off
% % semilogy(0:itout:maxit,DG_errorVec_CGNR,'-og','MarkerFaceColor','w','DisplayName','DG - CGN');
% % hold on;
% % semilogy(0:itout:maxit,DG_errorVec_GMRES,'-og','MarkerFaceColor','g','DisplayName','DG - GMRES');
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

% figure;
% loglog(hh,err,'-r');
% hold on
% loglog(hh,hh,'-b');
% loglog(hh,hh.^2,'-g');
% loglog(hh,hh.^3,'-y');
% loglog(hh,hh.^4,'-c');
% grid on
% legend('L^2-error', 'h', 'h^2', 'h^3', 'h^4');