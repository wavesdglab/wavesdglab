clear all;
%close all;

global k h

N=15;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);

computeSolNum2D = @computeSolNum2D_CG;

% Setup benchmark and parameters
benchmark = 'cavity';
switch benchmark
    case 'open'
        k = 15*pi;
        h = 1/16;
        tol = 1e-10; maxit = 1000; itout = 50;
    case 'cavity'
        k = 2.01*sqrt(2)*pi;
        h = 1/64;
        tol = 1e-10; maxit = 2000; itout =4;
        L = 1;
    case 'scatteringPML'
        k = 25;
        h = 0.05;
        tol = 1e-10; maxit = 2000; itout = 50;
        L = 1.1;
        R_disk = 1;
        L_PML = 0.2;
        computeSolNum2D = @computeSolNum2DPML_CG;
    case 'scattering_rect'
        k = 7.5*pi;
        h = 0.4;
        tol = 1e-7; maxit = 5000; itout = 10;
        L = 0.95;
        L_PML = 0.2;
        l = 0.5;
        computeSolNum2D = @computeSolNum2DPML_CG;
    case 'waveguide'
        k = 6*pi;
        h = 1/8;
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 1; % P1
PREC = 0; % for preconditioner

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree); % espace fonctionnel discret

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1); % nb de points par longueur d'onde

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CG - Benchmark "' benchmark '"']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

[~, sysA] = computeSolNum2D(mesh, dofm, PREC);

mat = sysA.matA;

% -------------------------------------------------------------------------
% Compute iterative solution
% -------------------------------------------------------------------------


[rrGMRES, ~, ~, ~, ~, X] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout,@computeNormError2D_CG);
solNum = X(:,end);


rrGMRES = [(0:itout:maxit)' rrGMRES];

errorL2 = computeNormError2D_CG(mesh, dofm, solNum);

solRef = computeSolProjL2_2D_CG(mesh, dofm);
errorProjL2 = computeNormError2D_CG(mesh, dofm, solRef);

disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

writeField2D(dofm, mesh, solNum, 'output/solNum.pos', "solNum");
writeField2D(dofm, mesh, solRef, 'output/solRef.pos', "solRef");
writeField2D(dofm, mesh, solNum-solRef, 'output/errNum.pos', "errNum");
system('gmsh benchmarks/cavity/cavity.msh output/solRef.pos output/solNum.pos output/errNum.pos&');



figure
hold on
set(0,'DefaultFigureWindowStyle','docked')

semilogy(0:itout:maxit,resVec,'b-o','DisplayName','Relative residual','linewidth', 1,'markersize', 5);
semilogy(0:itout:maxit,bound,'r--','DisplayName','Bound','linewidth', 1);
set(gca, 'YScale', 'log')
box on
grid on
xlim([0 maxit]);
ylim auto;
title(['CG - ' benchmark ' - GMRES - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)'], 'interpreter', 'latex', 'fontsize', 20)
xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
legend('Location', 'southwest', 'fontsize', 15)