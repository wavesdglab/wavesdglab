
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
tol = 1e-6; maxit = 2000; itout =4;
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


disp(['| Compute system and deflation subspace...']);

[~, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);

A = sysA.matA;
M = sysA.matP;
b = sysA.rhsA;
AMinv = A/M;

[eigvec,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,'closestEigvec',k);

[P,Q] = computeDefOp(nbEigVec, eigvec, A);

MinvP = M\P;

if maxit > size(A,2)
    maxit = size(A,2);
end


disp(['|             Done']);

% disp(['| Compute eigenvalues of A*(MinvP+Q)...']);
% [~, evdef] = eigs(A*(MinvP+Q),50-nbEigVec,'smallestabs');
% evdef = diag(evdef);
% evdef = [real(evdef) imag(evdef)];
% disp(['|             Done']);


% Compute GMRES with ADEF1 and closest eigvec and prec : A*(M\P+Q)*u = b, x = (M\P+Q)*u
disp(['| Preconditoned  GMRES with ADEF1...']);
sysA.matP = speye(size(A,1));
sysA.matA = A*(MinvP+Q);
sysA.rhsA = b;
[rrPAD, ~, itPAD, ~, ~, uPAD,hrv] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
xPAD = (MinvP+Q)*uPAD;
disp(['|             converges in ' num2str(itPAD) ' iterations']);



%%% Save results %%%

folder = "output/mainHistoryRitzValues_cavity";
if ~exist(folder, 'dir')
    mkdir(folder);
end

% for i=1:nbEigVec
%     writeField2D(dofm, mesh, eigvec(:,i), folder+"/eigvec"+num2str(i)+".pos", "eigvec"+num2str(i));
% end

% writeField2D(dofm, mesh, xPAD, folder+"/solADefP.pos", "solADefP");

csvwrite([folder+"/rrPAD.csv"],rrPAD);

csvwrite([folder+"/evdef.csv"],evdef);

csvwrite([folder+"/hrv.csv"],hrv);


%%% Plot results %%%

green = [0.4660 0.6740 0.1880];
magenta = [0.4940 0.1840 0.5560];
orange = [0.9290 0.6940 0.1250];
cyan = [0.3010 0.7450 0.9330];

%%% Spectra


% figure;
% hold on;

% s1 = scatter(real(evdef),imag(evdef),100, 'DisplayName','Eigenvalues of A*(M/P+Q)','linewidth', 2);
% s1.Marker = 'o';
% s1.MarkerEdgeColor = 'r';

% legend('Location', 'southwest', 'fontsize', 15)

% grid on; box on;
% title(['Spectra of $A$ - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)], 'interpreter', 'latex', 'fontsize', 20)
% set(0,'DefaultFigureWindowStyle','docked')



%%% Residuals

maxIt = max(itPAD);
minIt = min(itPAD);

disp(['ADEF1 and Shift : ' num2str(itPAD)]);

iterADefP = 0:itout:itout*size(rrPAD,1)-1;


figure
hold on
set(0,'DefaultFigureWindowStyle','docked')


p1  = semilogy(iterADefP,rrPAD,'y-o','DisplayName','Relative residual with ADEF1 and shift','linewidth', 2,'markersize', 10);
p1.Color = orange;

set(gca, 'YScale', 'log')
box on
grid on
xlim([0 maxIt+1]);
ylim auto;
title(['CG - ' benchmark ' - GMRES - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree) ' - nbEigvec=' num2str(nbEigVec)], 'interpreter', 'latex', 'fontsize', 20)
xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
legend('Location', 'southwest', 'fontsize', 15)