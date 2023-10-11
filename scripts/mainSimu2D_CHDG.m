clear all;
%close all;

% global k            % if the medium is homogeneous
global omega      % if the medium is heterogeneous

% Setup benchmark and parameters
benchmark = 'cavity_heterogeneous';
switch benchmark
    case 'open'
        k = 15*pi; %15*pi;
        h = 1/16;  % 1/16
        tol = 1e-10; maxit = 200; itout = 50;   %maxit = 1000
    case 'cavity'
        k = 7.1*sqrt(2)*pi; %7.01*sqrt(2)*pi;
        h = 0.1; %1/10;
        tol = 1e-10; maxit = 2000; itout = 100;
    case 'waveguide'
        k = 6*pi; %6*pi
        h = 1/6; %1/8
        tol = 1e-10; maxit = 4000; itout = 200;
    case 'cavity_heterogeneous'
        omega = 7.1*sqrt(2)*pi; %7.01*sqrt(2)*pi, 7.1*sqrt(2)*pi;
        h = 1/10; %1/10;
        tol = 1e-10; maxit = 2000; itout = 100;
end
degree = 5;
tau = 1;
BASIS = 0;
PREC = 1;
order=[1,1];  
% [1,1] or [1,2] first order transmission conditions and characteristic variables: standard CHDG
% [1,2]          first order transmission conditions and second order characteristic variables
% [2,2]          second order transmission conditions and characteristic variables

% Build mesh and DOF manager
mesh = benchmark2D(benchmark,h);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

%Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CHDG']);
disp(['---------------------------------------------------------']);
% disp(['    k                   ' num2str(k)]);
% disp(['    omega               ' num2str(omega)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    tau                 ' num2str(tau)]);
disp(['---------------------------------------------------------']);

% Compute numerical solution
% [solA, sysA] = computeSolNum2D_CHDG_new(mesh, dofm, tau, BASIS, PREC, order);
[solA, sysA] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, tau, BASIS, PREC);

% Compute numerical error
% errorL2 = computeNormError2D_DG(mesh, dofm, solA);

% Compute projection solution/errormesh.numTri*3*dofm.numDofPerTRI
% solP = computeSolProjL2_2D_DG(mesh, dofm);
% errorProjL2 = computeNormError2D_DG(mesh, dofm, solP);

% disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
% disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
% disp('---------------------------------------------------------');

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

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

% writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
% writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
% writeField2D(dofm, mesh, solA_1(1:mesh.numTri*3*dofm.numDofPerTRI)-solP, 'output/errNum.pos', "errNum");
% system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');
% system('gmsh output/solNum.pos&');

% writeField2D(dofm, mesh, solA_2, 'output/solNum.pos', "solNum");
% writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
% writeField2D(dofm, mesh, solA_2(1:mesh.numTri*3*dofm.numDofPerTRI)-solP, 'output/errNum.pos', "errNum");
% system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

% -------------------------------------------------------------------------
% Compute eigenvalues/eigenvectors
% -------------------------------------------------------------------------
% 
% mat = sysA.matPinv*sysA.matS;
% [eigenvec, eigenval] = eigs(mat,size(mat,1));
% eigenval = 1-diag(eigenval);
% rankEigenVec = rank(eigenvec);
% condEigenVec = cond(eigenvec);
% % 
% disp(['    Size                ' num2str(size(mat,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rankEigenVec)]);
% disp(['    Cond(eigenvectors)  ' num2str(condEigenVec)]);
% 
% figure;
% hold off; scatter(real(eigenval),imag(eigenval));
% axis equal
% %hold on; plot(fovals(mat,100));
% hold on; %plot(cos(0:0.01:2*pi)+1,sin(0:0.01:2*pi),'k');
% plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k');
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
        [resRedVec] = Richardson(sysA, tol, maxit, itout, alpha);
%         [resRedVec, resPhyVec, errorVec] = solverRichardson_DG(mesh, dofm, sysA, tol, maxit, itout, alpha, @computeNormError2D_DG);
    case 'CGNR'
        [resRedVec_1] = CGNR(sysA, tol, maxit, itout);
%         [resRedVec, resPhyVec, errorVec] = solverCGNRredu_DG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_DG);
    case 'GMRES'
        [resRedVec_1] = GMRES(sysA, tol, maxit, itout);
%         [resRedVec, resPhyVec, errorVec] = solverGMRESredu_DG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_DG);

end

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