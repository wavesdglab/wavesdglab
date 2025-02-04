
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
        h = 1/20;
        tol = 1e-10; maxit = 1000; itout = 50;
    case 'cavity'
%         k = 3.001*sqrt(2)*pi;
%           k = 3.1*sqrt(2)*pi;
        k = 5.877;
        h = 1/20;
        tol = 1e-8; maxit = 2000; itout =1;
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
        tol = 1e-7; maxit = 5000; itout = 1;
        L = 0.95;
        L_PML = 0.2;
        l = 0.5;
        computeSolNum2D = @computeSolNum2DPML_CG;
    case 'waveguide'
        k = 6*pi;
        h = 1/8;
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 2; % P1
PREC = 1; % for preconditioner

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

A = sysA.matA;
P = sysA.matP;
b = sysA.rhsA;

prec = 'right';
switch prec
    case 'left'
        [x1, ~, relres, it1, resVec1] = gmres(A,b,[],tol,maxit,eye(size(A,1)),P);
        [resVec2, ~, it2, ~, ~, x2] = solverGMRES_LP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
    case 'sym'
%         sqrtP = sqrtm(full(P));
        sysA.matP = eye(size(A,1))*9;
        sqrtP = eye(size(A,1))*3;
        [x1, ~, relres, it1, resVec1] = gmres(A,b,[],tol,maxit,sqrtP,sqrtP);
        [resVec2, ~, it2, ~, ~, x2] = solverGMRES_SP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
    case 'right'
%         [L,U] = ilu(P,struct('type','ilutp','droptol',1e-6));
        [L,U] = ilu(P);
%         [x1, ~, relres, it1, resVec1] = gmres(A/P,b,[],tol,maxit);
%         x1 = P\x1;
        [x1, ~, relres, it1, resVec1] = gmres(@(x) A*(U\(L\x)),b,[],tol,maxit);
        x1 = U\(L\x1);
        [resVec2, ~, it2, ~, ~, x2] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
    case 'none'
        sysA.matP = speye(size(A,1));
        [x1, ~, relres, it1, resVec1] = gmres(A,b,[],tol,maxit,eye(size(A,1)));
        [resVec2, ~, it2, ~, ~, x2] = solverGMRES_nP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
end

it1 = it1(2);
resVec1 = resVec1(:)./resVec1(1);

maxIt = min(it1,it2);

writeField2D(dofm, mesh, x1, 'output/x1.pos', "x1");
writeField2D(dofm, mesh, x2, 'output/x2.pos', "x2");
writeField2D(dofm, mesh, x1-x2, 'output/err.pos', "err");
system('gmsh output/mesh.msh output/x1.pos output/x2.pos output/err.pos&');

resVec2 = resVec2(1:maxIt+1);

error = norm(resVec1-resVec2);

disp([' ||r1-r2|| = ' num2str(error)]);


figure
hold on
set(0,'DefaultFigureWindowStyle','docked')

p1 = semilogy(0:itout:maxIt,resVec1,'b-o','DisplayName','Rel res (matlab)','linewidth', 2,'markersize', 5);
p2 = semilogy(0:itout:maxIt,resVec2,'r-x','DisplayName','Rel res (our gmres)','linewidth', 2,'markersize', 5);

set(gca, 'YScale', 'log')
box on
grid on
xlim([0 maxIt+1]);
ylim auto;
title(['CG - ' benchmark ' - GMRES ' prec ' preconditioning - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)], 'interpreter', 'latex', 'fontsize', 20)
xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
legend('Location', 'southwest', 'fontsize', 15)