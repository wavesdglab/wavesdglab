clear all;
%close all;

global k h

% Setup benchmark and parameters
benchmark = 'scattering_disk';
switch benchmark
    case 'open'
        k = 15*pi; h = 1/16;
        %k = 30*pi; h = 1/34;
        tol = 1e-6; maxit = 1000; itout = 50;
    case 'cavity'
        k = 7.1*sqrt(2)*pi; h = 1/10;
        %k = 7.01*sqrt(2)*pi; h = 1/15;
        tol = 1e-6; maxit = 2000; itout = 100;
    case 'waveguide'
        k = 6*pi; h = 1/8;
        %k = 12*pi; h = 1/17;
        tol = 1e-6; maxit = 4000; itout = 200;
    case 'scattering_disk'
        k = 2*pi; h = 1/10;
        %k = 30*pi; h = 1/34;
        tol = 1e-6; maxit = 1000; itout = 50;
end
tau = 1;
theta = 1;
degree = 3;
PREC = 0;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method DG']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['    tau                 ' num2str(tau)]);
disp(['    theta               ' num2str(theta)]);
disp(['---------------------------------------------------------']);

% Compute numerical solution/error
[solA, sysA] = computeSolNum2D_DG(mesh, dofm, tau, theta, PREC);
errorL2 = computeNormError2D_DG(mesh, dofm, solA);
disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);

% Compute projection solution/error
solP = computeSolProjL2_2D_DG(mesh, dofm);
errorProjL2 = computeNormError2D_DG(mesh, dofm, solP);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

% writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
% writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
% global PML_HIDE; PML_HIDE = 1;
% writeField2D(dofm, mesh, solA-solP, 'output/errNum.pos', "errNum");
% system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

% -------------------------------------------------------------------------
% Compute eigenvalues/eigenvectors
% -------------------------------------------------------------------------

% mat = sysA.matA;
% [eigenvec, eigenval] = eigs(mat,size(mat,1));
% eigenval = diag(eigenval);
% rankEigenVec = rank(eigenvec);
% condEigenVec = cond(eigenvec);
% 
% disp(['    Size                ' num2str(size(mat,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rankEigenVec)]);
% disp(['    Cond(eigenvectors)  ' num2str(condEigenVec)]);
% 
% % figure;
% % hold off; scatter(real(eigenval),imag(eigenval));
% % % hold on; plot(fovals(mat,100));
% % grid on; box on;
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
%     case 'CGNR'
%         [resVec, errorVec] = solverCGNR(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_DG);
%     case 'GMRES'
%         [resVec, errorVec] = solverGMRES(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_DG);
% end
% 
% figure;
% hold off
% semilogy(0:itout:maxit,resVec,'r','DisplayName','Relative residual');
% hold on
% semilogy(0:itout:maxit,errorVec,'k','DisplayName','Relative L2-error (iterative)');
% plot([0 maxit],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (direct)');
% plot([0 maxit],[errorProjL2 errorProjL2],'k:','DisplayName','Relative L2-error (projection)');
% box on;
% grid on;
% legend('Location','southwest');
% title(['DG - ' benchmark ' - ' solver ' - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)])
% xlim([0 maxit]);
% ylim([0.005 1]);
% xlabel('Iteration');
% ylabel('Value');