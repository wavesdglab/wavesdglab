
clear all;
%close all;

global k h

N=15;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);


% Setup benchmark and parameters
benchmark = 'cavity';
k = 3.01*sqrt(2)*pi;
h = 1/64;
maxit = 2000; itout =4;
L = 1;
degree = 1;
PREC = 0;
nbEigVec=1;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

% -------------------------------------------------------------------------
% Compute solution
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CG - Benchmark "' benchmark '"']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

[~, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);

tabTol = [1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10];

nbitGMRES = zeros(2, length(tabTol));

nbitAD = zeros(2, length(tabTol));
nbitD = zeros(2, length(tabTol));

nbitGMRES(1,:) = tabTol';
nbitAD(1,:) = tabTol';
nbitD(1,:) = tabTol';


[eigvec,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,'closestEigvec',k);

A = sysA.matA;
M = sysA.matP;
b = sysA.rhsA;
AMinv = A/M;

if maxit > size(A,2)
    maxit = size(A,2);
end

[P,Q] = computeDefOp(nbEigVec, eigvec, A);

MinvP = M\P;

for tol = tabTol
    
    
    %%%%%%%%%%% No deflation %%%%%%%%%%%
    
    % Compute GMRES without prec
    disp(['| GMRES...']);
    [xGMRES, ~, ~, itGMRES, rrGMRES] = gmres(A,b,[],tol,maxit);
    itGMRES = itGMRES(2);
    disp(['|             converges in ' num2str(itGMRES) ' iterations']);
    
    nbitGMRES(2, tabTol == tol) = itGMRES;
    
    % Compute GMRES with ADEF1 and closest eigvec : A*(P+Q)*u = b, x = (P+Q)*u
    disp(['| GMRES with ADEF1...']);
    [uAD, ~, ~, itAD, rrAD] = gmres(A*(P+Q),b,[],tol,maxit);
    itAD = itAD(2);
    disp(['|             converges in ' num2str(itAD) ' iterations']);
    
    nbitAD(2, tabTol == tol) = itAD;
    
    % Compute GMRES with DEF1 and closest eigvec : P*A*x = P*b
    disp(['| GMRES with DEF1...']);
    [xD, ~, ~, itD, rrD] = gmres(P*A,P*b,[],tol,maxit);
    itD = itD(2);
    disp(['|             converges in ' num2str(itD) ' iterations']);

    nbitD(2, tabTol == tol) = itD;
    
    
    disp(['---------------------------------------------------------']);
    disp(['Number of iterations GMRES: ' num2str(itGMRES)]);
    disp(['Number of iterations GMRES with ADEF: ' num2str(itAD)]);
    disp(['Number of iterations GMRES with DEF: ' num2str(itD)]);
    disp(['---------------------------------------------------------']);
    
    
    
    
end

csvwrite('output/nbitGMRES.csv', nbitGMRES');
csvwrite('output/nbitAD.csv', nbitAD');
csvwrite('output/nbitD.csv', nbitD');

% maxIt = max(max(nbitGMRES), max(max(nbitAD), max(nbitD)));
% minIt = min(min(nbitGMRES), min(min(nbitAD), min(nbitD)));

maxIt = max(max(nbitGMRES(2,:)), max(max(nbitAD(2,:)), max(nbitD(2,:))));
minIt = min(min(nbitGMRES(2,:)), min(min(nbitAD(2,:)), min(nbitD(2,:))));


figure
hold on
set(0,'DefaultFigureWindowStyle','docked')

p1 = semilogy(tabTol,nbitGMRES(2,:),'b-o','DisplayName','rrGMRES','linewidth', 2,'markersize', 10);
p2 = semilogy(tabTol,nbitAD(2,:),'r-+','DisplayName','rrAD','linewidth', 2,'markersize', 10);
p3 = semilogy(tabTol,nbitD(2,:),'g-x','DisplayName','rrD','linewidth', 2,'markersize', 10);

set(gca, 'YScale', 'log')
box on
grid on
ylim auto;
title(['CG - ' benchmark ' - GMRES - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree) ' - nbEigvec=' num2str(nbEigVec)], 'interpreter', 'latex', 'fontsize', 20)
xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
legend('Location', 'southwest', 'fontsize', 15)
