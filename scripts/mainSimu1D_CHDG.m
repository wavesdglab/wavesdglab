clear all;
%close all;

global k BCLeft BCRight

% Define parameters
degree = 3;
k = 40;
numE = 200;
tau = 1; % 1i
resTol = 1e-4;
PREC = 'PrecNone'; % PrecNone PrecMass PrecMass2 PrecDiag
BCLeft = 'DIR';
BCRight = 'ABC';

% Build mesh and dofManager
mesh = buildMesh1D(0, 1, numE);
dofm = buildDofManager1D_DG(mesh, degree);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CHDG - Benchmark ' BCLeft '/' BCRight ]);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(1/numE)]);
disp(['    degree              ' num2str(degree)]);
disp(['    numE                ' num2str(numE)]);
disp(['    tau                 ' num2str(tau)]);
disp(['---------------------------------------------------------']);

[solA, sysA] =  computeSolNum1D_CHDG(mesh, dofm, tau);
%[solA, sysA] =  computeSolNum1D_CHDGb(mesh, dofm, tau);
errorL2 = computeNormError1D_DG(mesh, dofm, solA);

solP = computeSolProjL2_1D_DG(mesh, dofm);
errorProjL2 = computeNormError1D_DG(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.6e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.6e')]);
disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Vizu solution
% -------------------------------------------------------------------------

plotField1D(mesh, dofm, solA, 'Numerical solution');

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
% % Plot spectrum
% figure;
% hold off; scatter(real(eigenval),imag(eigenval));
% % hold on; plot(fovals(mat,100));
% hold on; plot(cos(0:0.01:2*pi)+1,sin(0:0.01:2*pi),'k');
% grid on; box on;
% 
% % Compute condition number
% condestMat = condest(mat);
% 
% disp(['    Cond(mat)           ' num2str(condestMat,'%1.2e')]);
% disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Compute iterative solution
% -------------------------------------------------------------------------

solver = 'GMRES';
switch solver
    case 'CGN'
        tol = 1e-10; maxit = 1000; itout = 10;
        [resRedVec, resPhyVec, errorVec, iter, flag] = solverCGNredu_DG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError1D_DG);
    case 'GMRES'
        tol = 1e-10; maxit = numE; itout = 1;
        [resRedVec, resPhyVec, errorVec, iter, flag] = solverGMRESredu_DG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError1D_DG);
end

figure;
hold off
semilogy(0:itout:maxit,resPhyVec,'r-o','DisplayName','Relative residual (Phy)');
hold on
semilogy(0:itout:maxit,resRedVec,'b-','DisplayName','Relative residual (Red)');
semilogy(0:itout:maxit,errorVec,'k','DisplayName','Relative L2-error (iterative)');
plot([0 maxit],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (direct)');
plot([0 maxit],[errorProjL2 errorProjL2],'k:','DisplayName','Relative L2-error (projection)');
box on;
grid on;
legend('Location','southwest');
title(['CG - ' BCLeft '/' BCRight ' - ' solver ' - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)])
xlim([0 maxit]);
xlabel('Iteration');
ylabel('Value');
