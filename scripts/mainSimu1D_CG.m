clear all;
%close all;

global k BCLeft BCRight

% Setup benchmark and parameters
degree = 3;
k = 4;
numE = 100;
h = 1/numE;
resTol = 1e-4;
PREC = 'PrecNone'; % PrecNone PrecMass PrecDiag PrecShiftLap
shiftPrec = 1; % 1+1i
BCLeft = 'ABC';
BCRight = 'ABC';

% Build mesh and DOF manager
mesh = buildMesh1D(0, 1, numE);
dofm = buildDofManager1D_CG(mesh, degree);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CG - Benchmark ' BCLeft '/' BCRight ]);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    numE                ' num2str(numE)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum1D_CG(mesh, dofm, PREC, shiftPrec);
errorL2 = computeNormError1D_CG(mesh, dofm, solA);

solP = computeSolProjL2_1D_CG(mesh, dofm);
errorProjL2 = computeNormError1D_CG(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.6e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.6e')]);
disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Vizu solution
% -------------------------------------------------------------------------

% figure;
% plotField1D(mesh, dofm, solA, 'Numerical solution');

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
        [resRedVec, resPhyVec, errorVec, iter, flag] = solverCGNredu_CG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError1D_CG);
    case 'GMRES'
        tol = 1e-10; maxit = numE; itout = 1;
        [resRedVec, resPhyVec, errorVec, iter, flag] = solverGMRESredu_CG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError1D_CG);
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
