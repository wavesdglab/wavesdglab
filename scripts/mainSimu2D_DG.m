%close all;
clear all;

global k

% Setup benchmark and parameters
benchmark = 'open';
switch benchmark
    case 'open'
        k = 15*pi;
        h = 1/16;
    case 'cavity'
        k = 7.1*sqrt(2)*pi;
        h = 1/10;
    case 'waveguide'
        k = 6*pi;
        h = 1/8;
end
degree = 3;
tau = 1;
theta = 0;

% Build mesh and DOF manager
mesh = benchmark2D(benchmark,h);
mesh = buildMeshConnectivity(mesh);
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
[solA, sysA] = computeSolNum2D_DG(mesh, dofm, tau, theta);
errorL2 = computeNormError2D_DG(mesh, dofm, solA);

% Compute projection solution/error
solP = computeSolProjL2_2D_DG(mesh, dofm);
errorProjL2 = computeNormError2D_DG(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
disp(['---------------------------------------------------------']);

% disp([num2str(k) ' ' num2str(h) ' ' num2str(degree) ' ' num2str(Dlambda) ' ' num2str(errorL2) ' ' num2str(errorProjL2)]);

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

% writeField(dofm, mesh, solA, 'output/solNum.pos', "solNum");
% writeField(dofm, mesh, solP, 'output/solRef.pos', "solRef");
% writeField(dofm, mesh, solA-solP, 'output/errNum.pos', "errNum");
% system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

% -------------------------------------------------------------------------
% Compute eigenvalues/eigenvectors
% -------------------------------------------------------------------------

% mat = sysA.matA;
% [eigenvec,eigenval] = eigs(mat,size(mat,1));
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

solver = 'CGN';
switch solver
    case 'CGN'
        tol = 1e-10; maxit = 1000; itout = 10;
        [resVec, errorVec, iter, flag] = solverCGN_DG(mesh, dofm, sysA, tol, maxit, itout);
    case 'GMRES'
        tol = 1e-10; maxit = 50; itout = 1;
        [resVec, errorVec, iter, flag] = solverGMRES_DG(mesh, dofm, sysA, tol, maxit, itout);
end

figure;
hold off
semilogy(0:itout:maxit,resVec,'r','DisplayName','Relative residual');
hold on
semilogy(0:itout:maxit,errorVec,'k','DisplayName','Relative L2-error (iterative)');
plot([0 maxit],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (direct)');
plot([0 maxit],[errorProjL2 errorProjL2],'k:','DisplayName','Relative L2-error (projection)');
box on;
grid on;
legend('Location','southwest');
title(['DG - ' benchmark ' - ' solver ' - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)])
xlim([0 maxit]);
ylim([0.005 1]);
xlabel('Iteration');
ylabel('Value');
