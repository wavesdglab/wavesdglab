
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
PREC = 1;
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
    M = sysA.matP;
    
    if maxit > size(A,2)
        maxit = size(A,2);
    end

    [Pdef,Qdef,Q] = computeDefOp(nbEigVec, eigvec, A);

    AMinv = A/M;

    
    %%%%%%%%%%% No deflation %%%%%%%%%%%
    
    % Compute GMRES without prec
    [xGMRES, ~, ~, itGMRES, rrGMRES] = gmres(A,b,[],tol,maxit);
    itGMRES = itGMRES(2);
    nbit(tabfreq == k, 2) = itGMRES;
    
    
    % Compute GMRES with right preconditioner : A*Minv*u = b, x = Minv*u
    [xGMRESPRight, ~, ~, itGMRESPRight, rrGMRESPRight] = gmres(AMinv,b,[],tol,maxit);
    itGMRESPRight = itGMRESPRight(2);
    nbit(tabfreq == k, 3) = itGMRESPRight;
    
    
    % Compute GMRES with left Def : Pdef*A*u = Pdef*b, x = Q*b + Qdef*u
    [uDLeft, ~, ~, itDLeft, rrDLeft] = gmres(Pdef*A,Pdef*b,[],tol,maxit);
    itDLeft = itDLeft(2);
    nbit(tabfreq == k, 4) = itDLeft;


    % Compute GMRES with left Def and right prec : Pdef*AMinv*u = Pdef*b, x = Q*b + Qdef*Minv*u
    [uDLeftPRight, ~, ~, itDLeftPRight, rrDLeftPRight] = gmres(Pdef*AMinv,Pdef*b,[],tol,maxit);
    itDLeftPRight = itDLeftPRight(2);
    nbit(tabfreq == k, 5) = itDLeftPRight;


    % Compute GMRES with additive right method : A*(I+Q)*u = b, x = (I+Q)*u
    [xAddRight, ~, ~, itAddRight, rrAddRight] = gmres(A*(eye(size(A))+Q),b,[],tol,maxit);
    itAddRight = itAddRight(2);
    nbit(tabfreq == k, 6) = itAddRight;


    % Compute GMRES with additive right method and right preconditioner : A*(Minv+Q)*u = b, x = (Minv+Q)*u
    [xAddRightPrec, ~, ~, itAddRightPrec, rrAddRightPrec] = gmres(AMinv+A*Q,b,[],tol,maxit);
    itAddRightPrec = itAddRightPrec(2);
    nbit(tabfreq == k, 7) = itAddRightPrec;

    
    disp(['---------------------------------------------------------']);
    disp(['Number of iterations GMRES: ' num2str(itGMRES)]);
    disp(['Number of iterations GMRES with right preconditioner: ' num2str(itGMRESPRight)]);
    disp(['Number of iterations GMRES with left deflation: ' num2str(itDLeft)]);
    disp(['Number of iterations GMRES with left deflation and right preconditioner: ' num2str(itDLeftPRight)]);
    disp(['Number of iterations GMRES with additive right method: ' num2str(itAddRight)]);
    disp(['Number of iterations GMRES with additive right method and right preconditioner: ' num2str(itAddRightPrec)]);
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

