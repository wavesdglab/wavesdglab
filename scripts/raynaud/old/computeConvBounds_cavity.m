
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

k = 3.1*sqrt(2)*pi;
h = 1/32;
tol = 1e-8; maxit = 2000; itout =4;
L = 1;

degree = 1;
PREC = 1;
nbEigVec=30;

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

A = sysA.matA;
M = sysA.matP;
b = sysA.rhsA;

if maxit > size(A,2)
    maxit = size(A,2);
end


[~, evA] = eigs(A,size(A,1),'smallestabs');
evA = diag(evA);

[~, evMA] = eigs(M\A,size(A,1),'smallestabs');
evMA = diag(evMA);



%%%%%%%%%%% No deflation %%%%%%%%%%%

% Compute GMRES with prec
[xGMRESP, ~, ~, itGMRESP, rrGMRESP] = gmres(A/M,b,[],tol,maxit);
itGMRESP = itGMRESP(2);
rrGMRESP = rrGMRESP(:)./rrGMRESP(1);
rrGMRESP = rrGMRESP(1:itout:end);
xGMRESP = M\xGMRESP;



% Compute GMRES without prec
[xGMRES, ~, ~, itGMRES, rrGMRES] = gmres(A,b,[],tol,maxit);
itGMRES = itGMRES(2);
rrGMRES = rrGMRES(:)./rrGMRES(1);
rrGMRES = rrGMRES(1:itout:end);



%%%%%%%%%%% Deflation %%%%%%%%%%%




[~,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,"firstEigvec");
[eigvec,~] = eigs(A,nbEigVec,'smallestabs');

[P,Q] = computeDefOp(nbEigVec, eigvec, A);

[~, evdef] = eigs((P+Q)*A,size(A,1),'smallestabs');
evdef = diag(evdef);



%%% No preconditioner :

% Compute GMRES with DEF1 and closest eigvec : P*A*x = P*b
[xD, ~, ~, itD, rrD] = gmres(P*A,P*b,[],tol,maxit);
itD = itD(2);
rrD = rrD(:)./rrD(1);
rrD = rrD(1:itout:end);
xD = P'*xD + Q*b;



% Compute GMRES with ADEF1 and closest eigvec : (P+Q)*A*x = (P+Q)*b
[xAD, ~, ~, itAD, rrAD] = gmres((P+Q)*A,(P+Q)*b,[],tol,maxit);
itAD = itAD(2);
rrAD = rrAD(:)./rrAD(1);
rrAD = rrAD(1:itout:end);
xAD = (P+Q)*xAD;



%%% Add preconditioner :

% Compute GMRES with DEF1 and closest eigvec and prec : P*A/M*u = P*b, x = M\u
[xPD, ~, ~, itPD, rrPD] = gmres(P*A/M,P*b,[],tol,maxit);
itPD = itPD(2);
rrPD = rrPD(:)./rrPD(1);
rrPD = rrPD(1:itout:end);
xPD = P'/M*xPD + Q*b;



% Compute GMRES with ADEF1 and closest eigvec and prec : A*(M\P+Q)u = b, x = (M\P+Q)u
MinvP_Q = M\P+Q;
[xPAD, ~, ~, itPAD, rrPAD] = gmres(A*MinvP_Q,b,[],tol,maxit);
itPAD = itPAD(2);
rrPAD = rrPAD(:)./rrPAD(1);
rrPAD = rrPAD(1:itout:end);
xPAD = MinvP_Q*xPAD;


%%%%%%%%%%% Computing convergence bounds %%%%%%%%%%%


lambda_min = min(abs(evdef));
lambda_max = max(abs(evdef));
kappa = lambda_max/lambda_min;

% disp(['Condition number of A*(P+Q) : ' num2str(kappa)]);
rho = (lambda_max-lambda_min)/2;
mu = rho + lambda_min;

bound = zeros(1,size(rrAD,1));
for i = 1:size(rrAD,1)
    bound(i) = 2*(kappa-1)/(kappa+1)^(i/2);
%     bound(i) = (rho/mu)^(i);
end

%%% Plot results %%%

green = [0.4660 0.6740 0.1880];
magenta = [0.4940 0.1840 0.5560];
orange = [0.9290 0.6940 0.1250];
cyan = [0.3010 0.7450 0.9330];

%%% Spectra


figure;
hold on;
s1 = scatter(real(evA),imag(evA),100,'DisplayName','Eigenvalues of A');
s1.Marker = '+';
s1.MarkerEdgeColor = 'b';

s2 = scatter(real(evdef),imag(evdef),100, 'DisplayName','Eigenvalues of (P+Q)*A');
s2.Marker = 'o';
s2.MarkerEdgeColor = 'r';

legend('Location', 'southwest', 'fontsize', 15)

grid on; box on;
title(['Spectra of $A$ - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)], 'interpreter', 'latex', 'fontsize', 20)
% xlim([-0.04 0.06]);
% hold on; plot(fovals(A));
set(0,'DefaultFigureWindowStyle','docked')



%%% Residuals

it = [itGMRES itGMRESP itD itAD itPD itPAD];

maxIt = max(it);
minIt = min(it);

disp(['GMRES: ' num2str(itGMRES)]);
disp(['GMRES and Shift: ' num2str(itGMRESP)]);
disp(['DEF1 : ' num2str(itD)]);
disp(['ADEF1 : ' num2str(itAD)]);
disp(['DEF1 and Shift : ' num2str(itPD)]);
disp(['ADEF1 and Shift : ' num2str(itPAD)]);

disp(['Difference: ' num2str(100*(maxIt-minIt)/maxIt) '%']);


iterGMRES = 0:itout:itout*size(rrGMRES,1)-1;
iterGMRESP = 0:itout:itout*size(rrGMRESP,1)-1;
iterDef = 0:itout:itout*size(rrD,1)-1;
iterADef = 0:itout:itout*size(rrAD,1)-1;
iterDefP = 0:itout:itout*size(rrPD,1)-1;
iterADefP = 0:itout:itout*size(rrPAD,1)-1;


figure
hold on
set(0,'DefaultFigureWindowStyle','docked')

p1 = semilogy(iterGMRES,rrGMRES,'b-o','DisplayName','Relative residual','linewidth', 2,'markersize', 10);
p2 = semilogy(iterGMRESP,rrGMRESP,'r-o','DisplayName','Relative residual with shift','linewidth', 2,'markersize', 10);
p3 = semilogy(iterDef,rrD,'g-o','DisplayName','Relative residual with DEF1','linewidth', 2,'markersize', 10);
p3.Color = green;
p4 = semilogy(iterADef,rrAD,'c-o','DisplayName','Relative residual with ADEF1','linewidth', 2,'markersize', 10);
p4.Color = cyan;
p5 = semilogy(iterDefP,rrPD,'m-o','DisplayName','Relative residual with DEF1 and shift','linewidth', 2,'markersize', 10);
p5.Color = magenta;
p6  = semilogy(iterADefP,rrPAD,'y-o','DisplayName','Relative residual with ADEF1 and shift','linewidth', 2,'markersize', 10);
p6.Color = orange;

p7 = semilogy(iterADef,bound,'--','DisplayName','Bound','linewidth', 2);

set(gca, 'YScale', 'log')
box on
grid on
xlim([0 maxIt+1]);
ylim auto;
title(['CG - ' benchmark ' - GMRES - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree) ' - nbEigvec=' num2str(nbEigVec)], 'interpreter', 'latex', 'fontsize', 20)
xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
legend('Location', 'southwest', 'fontsize', 15)