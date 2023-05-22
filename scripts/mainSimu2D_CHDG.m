clear all;
%close all;
tic
global k

% for i = 0:5

% Setup benchmark and parameters
benchmark = 'cavity';
switch benchmark
    case 'open'
        k = 15*pi; %15*pi;
        h = 1/16;  % 1/16
%         h = 1/5*(1/2)^i;
%         H(i+1)=h;
        tol = 1e-10; maxit = 1000; itout = 50;
    case 'cavity'
        k = 7.1*sqrt(2)*pi; %7.01*sqrt(2)*pi;
        h = 1/10; %1/10;
%         h = 1/5*(1/2)^i;
%         H(i+1)=h;
        tol = 1e-10; maxit = 2000; itout = 100;
    case 'waveguide'
        k = 6*pi; %6*pi
        h = 1/8; %1/8
%         h = 1/5*(1/2)^i;
%         H(i+1)=h;
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 3;
tau = 1;
BASIS = 0;
PREC = 1; 

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

% Compute numerical solution

[solA_1, sysA_1] = computeSolNum2D_CHDG_1st_order(mesh, dofm, tau, BASIS, PREC);
% [solA_1, sysA_1] = computeSolNum2D_CHDG_1st_order_ALT(mesh, dofm, tau, BASIS, PREC);

% [solA_2, sysA_2] = computeSolNum2D_CHDG_2nd_order_ALT(mesh, dofm, tau, BASIS, PREC);
% [solA_2, sysA_2] = computeSolNum2D_CHDG_2nd_order(mesh, dofm, tau, BASIS, PREC);

% CHDG2
[solA_2, sysA_2] = computeSolNum2D_CHDG2_2nd_order_FULL(mesh, dofm, tau, BASIS, PREC);
% [solA_2, sysA_2] = computeSolNum2D_CHDG2_2nd_order(mesh, dofm, tau, BASIS, PREC);

% CHDG3
% [solA_2, sysA_2] = computeSolNum2D_CHDG3_2nd_order_FULL(mesh, dofm, tau, BASIS, PREC);
% [solA_2, sysA_2] = computeSolNum2D_CHDG3_2nd_order(mesh, dofm, tau, BASIS, PREC);                     % IT WORKS

% Compute numerical error

errorL2_1 = computeNormError2D_DG(mesh, dofm, solA_1);
errorL2_2 = computeNormError2D_DG(mesh, dofm, solA_2);

% Compute projection solution/errormesh.numTri*3*dofm.numDofPerTRI
solP = computeSolProjL2_2D_DG(mesh, dofm);
errorProjL2 = computeNormError2D_DG(mesh, dofm, solP);
% % % 
disp('     1st order');
disp(['    L2-Error (numSol)   ' num2str(errorL2_1,'%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
disp('---------------------------------------------------------');
% 
% err(i+1,1) = errorL2;
% err(i+1,3) = errorProjL2;

% [solA, sysA] = computeSolNum2D_CHDG_2nd_order(mesh, dofm, tau, BASIS, PREC);
% errorL2 = computeNormError2D_DG(mesh, dofm, solA);
% 
disp('     2nd order');
disp(['    L2-Error (numSol)   ' num2str(errorL2_2,'%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
disp('---------------------------------------------------------');
% 
% err(i+1,2) = errorL2;
% 
% end
% 
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

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

% writeField2D(dofm, mesh, solA_2, 'output/solNum.pos', "solNum");
% writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
% writeField2D(dofm, mesh, solA_2(1:mesh.numTri*3*dofm.numDofPerTRI)-solP, 'output/errNum.pos', "errNum");
% system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

% -------------------------------------------------------------------------
% Compute eigenvalues/eigenvectors
% -------------------------------------------------------------------------

% mat1 = sysA_1.matPinv*sysA_1.matS;
% [eigenvec1, eigenval1] = eigs(mat1,size(mat1,1));
% eigenval1 = 1-eigenval1;
% eigenval1 = diag(eigenval1);
% rankEigenVec1 = rank(eigenvec1);
% condEigenVec1 = cond(eigenvec1);

% disp(['    Size                ' num2str(size(mat1,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rankEigenVec1)]);
% disp(['    Cond(eigenvectors)  ' num2str(condEigenVec1)]);
% 
% figure;
% hold off; scatter(real(eigenval1),imag(eigenval1));
% axis equal
% %hold on; plot(fovals(mat,100));
% hold on; %plot(cos(0:0.01:2*pi)+1,sin(0:0.01:2*pi),'k');
% plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k');
% grid on; box on;
% set(gcf, 'PaperUnits', 'points','PaperPosition', [0 0 500 500]);
% print(['output/Eigenvalues-' benchmark '-CHDG.png'],'-dpng');

% mat2 = sysA_2.matPinv*sysA_2.matS;
% [eigenvec2, eigenval2] = eigs(mat2,size(mat2,1));
% eigenval2 = diag(eigenval2);
% rankEigenVec2 = rank(eigenvec2);
% condEigenVec2 = cond(eigenvec2);
% 
% disp(['    Size                ' num2str(size(mat2,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rankEigenVec2)]);
% disp(['    Cond(eigenvectors)  ' num2str(condEigenVec2)]);
% 
% figure;
% hold off; scatter(real(eigenval2),imag(eigenval2));
% axis equal
% %hold on; plot(fovals(mat,100));
% hold on; plot(cos(0:0.01:2*pi)+1,sin(0:0.01:2*pi),'k');
% grid on; box on;
% set(gcf, 'PaperUnits', 'points','PaperPosition', [0 0 500 500]);
% print(['output/Eigenvalues-' benchmark '-CHDG.png'],'-dpng');

% % Compute condition number
% condestMat = condest(mat);
% 
% disp(['    Cond(mat)           ' num2str(condestMat,'%1.2e')]);
% disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Compute iterative solution
% -------------------------------------------------------------------------

solver = 'Rich';
switch solver
    case 'Rich'
        alpha = 1;
        [resRedVec_1, resPhyVec_1, errorVec_1] = solverRichardson_DG(mesh, dofm, sysA_1, tol, maxit, itout, alpha, @computeNormError2D_DG);
        [resRedVec_2, resPhyVec_2, errorVec_2] = solverRichardson_DG(mesh, dofm, sysA_2, tol, maxit, itout, alpha, @computeNormError2D_DG);
    case 'CGNR'
        [resRedVec_1, resPhyVec_1, errorVec_1] = solverCGNRredu_DG(mesh, dofm, sysA_1, tol, maxit, itout, @computeNormError2D_DG);
        [resRedVec_2, resPhyVec_2, errorVec_2] = solverCGNRredu_DG(mesh, dofm, sysA_2, tol, maxit, itout, @computeNormError2D_DG);
    case 'GMRES'
        [resRedVec_1, resPhyVec_1, errorVec_1] = solverGMRESredu_DG(mesh, dofm, sysA_1, tol, maxit, itout, @computeNormError2D_DG);
        [resRedVec_2, resPhyVec_2, errorVec_2] = solverGMRESredu_DG(mesh, dofm, sysA_2, tol, maxit, itout, @computeNormError2D_DG);
end

% errorVec(2*j+1,:) = errorVec_1;
% errorVec(2*j+2,:) = errorVec_2;

errorVec(3,:) = errorVec_1;
errorVec(4,:) = errorVec_2;

% end

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

figure(9);
hold off
semilogy(0:itout:maxit,errorVec(3,:),'r-x','DisplayName','Relative L2-error (iterative - 1st order) - precond.');
hold on
semilogy(0:itout:maxit,errorVec(4,:),'g-x','DisplayName','Relative L2-error (iterative - 2nd order) - precond.');
hold on
plot([0 maxit],[errorL2_1 errorL2_1],'k--','DisplayName','Relative L2-error (direct - 1st order)');
hold on
plot([0 maxit],[errorL2_2 errorL2_2],'k-.','DisplayName','Relative L2-error (direct - 2nd order)');
box on;
grid on;
% legend('Location','southwest');
legend('Iterative - 1st order - precond.', 'Iterative - 2nd order - precond.', 'Direct - 1st order', 'Direct - 2nd order');
title(['CHDG - ' benchmark ' - ' solver ' - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)])
xlim([0 maxit]);
ylim([0.005 1]);
xlabel('Iteration');
ylabel('Relative L^2-error');

toc