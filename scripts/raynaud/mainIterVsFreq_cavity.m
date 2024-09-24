
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
h = 1/64;
tol = 1e-6; maxit = 2000; itout =4;
L = 1;
degree = 1;
PREC = 0;
nbEigVec=1;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

tabfreq = 2.851*sqrt(2)*pi:0.001*sqrt(2)*pi:3.251*sqrt(2)*pi;


nbit = zeros(length(tabfreq), 4);
nbit(:,1) = tabfreq;

[eigvec,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,'closestEigvec',3.00001*sqrt(2)*pi);

for k = tabfreq
    
    
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
    
    A = sysA.matA;
    b = sysA.rhsA;
    
    if maxit > size(A,2)
        maxit = size(A,2);
    end

    [P,Q] = computeDefOp(nbEigVec, eigvec, A);
    
    %%%%%%%%%%% No deflation %%%%%%%%%%%
    
    % Compute GMRES without prec
    [xGMRES, ~, ~, itGMRES, rrGMRES] = gmres(A,b,[],tol,maxit);
    itGMRES = itGMRES(2);
    nbit(tabfreq == k, 2) = itGMRES;
    
    
    % Compute GMRES with ADEF1 and closest eigvec : (P+Q)*A*x = (P+Q)*b
    [uAD, ~, ~, itAD, rrAD] = gmres(A*(P+Q),b,[],tol,maxit);
    itAD = itAD(2);
    nbit(tabfreq == k, 3) = itAD;
    
    
    % Compute GMRES with DEF1 and closest eigvec : P*A*x = P*b
    [xD, ~, ~, itD, rrD] = gmres(P*A,P*b,[],tol,maxit);
    itD = itD(2);
    nbit(tabfreq == k, 4) = itD;
    
    disp(['---------------------------------------------------------']);
    disp(['Number of iterations GMRES: ' num2str(itGMRES)]);
    disp(['Number of iterations GMRES with ADEF: ' num2str(itAD)]);
    disp(['Number of iterations GMRES with DEF: ' num2str(itD)]);
    disp(['---------------------------------------------------------']);
    
    
    
    
end

csvwrite('output/nbit.csv',nbit);

% maxIt = max(nbit(:,2:4));
% minIt = min(nbit(:,2:4));

% figure
% hold on
% set(0,'DefaultFigureWindowStyle','docked')

% p1 = semilogy(tabfreq,nbit(:,2),'b-o','DisplayName','GMRES','linewidth', 2,'markersize', 10);
% p2 = semilogy(tabfreq,nbit(:,3),'r-+','DisplayName','rrAD','linewidth', 2,'markersize', 10);
% p3 = semilogy(tabfreq,nbit(:,4),'g-x','DisplayName','rrD','linewidth', 2,'markersize', 10);

% set(gca, 'YScale', 'log')
% box on
% grid on
% xlim([2.5*sqrt(2)*pi 3.5*sqrt(2)*pi]);
% ylim auto;
% title(['CG - ' benchmark ' - GMRES - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree) ' - nbEigvec=' num2str(nbEigVec)], 'interpreter', 'latex', 'fontsize', 20)
% xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
% ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
% legend('Location', 'southwest', 'fontsize', 15)

