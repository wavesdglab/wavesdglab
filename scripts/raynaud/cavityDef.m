
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
        k = 3.01*sqrt(2)*pi;
        h = 1/64;
        tol = 1e-6; maxit = 2000; itout =4;
        L = 1;
    case 'scatteringPML'
        k = 25;
        h = 0.05;
        tol = 1e-10; maxit = 2000; itout = 50;
        L = 1.1;
        R_disk = 1;
        L_PML = 0.2;
        computeSolNum2D = @computeSolNum2DPML_CG;
    case 'scattering_rec'
        global LdomX LdomY LpmlX LpmlY
        k = 10.5*pi;
        h = 1/8;
        tol = 1e-7; maxit = 5000; itout = 5;
        LdomX = 0.95;
        LdomY = 0.5;
        LpmlX = 0.8;
        LpmlY = 0.8;
    case 'waveguide'
        k = 6*pi;
        h = 1/8;
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 1; % P1
PREC = 1; % for preconditioner
eigvecToDeflate = "closesteigvec"; %"firsteigvec" or "closesteigvec"
nbEigVec=5;

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

disp(['| Compute system...']);

[~, sysA] = computeSolNum2D(mesh, dofm, PREC);

disp(['|             Done']);

A = sysA.matA;
M = sysA.matP;
b = sysA.rhsA;

if maxit > size(A,2)
    maxit = size(A,2);
end

% disp(['| Compute eigenvalues of A...']);
% 
% [~, evA] = eigs(A,size(A,2),'smallestabs');
% evA = diag(evA);
% evA = [real(evA) imag(evA)];
% csvwrite(["output/evA.csv"],evA);

% disp(['|             Done']);

% disp(['| Compute A/M...']);

AMinv = A/M;
% disp(['|             Done']);
% disp(['| Compute eigenvalues of A/M...']);
% [~, evAMinv] = eigs(AMinv,50,'smallestabs');
% evAMinv = diag(evAMinv);
% evAMinv = [real(evAMinv) imag(evAMinv)];
% csvwrite(["output/spectrum_prec.csv"],evAMinv);
% disp(['|             Done']);


%%%%%%%%%%% No deflation %%%%%%%%%%%

% Compute GMRES with prec
disp(['| Preconditioned GMRES...']);
[xGMRESP, ~, ~, itGMRESP, rrGMRESP] = gmres(AMinv,b,[],tol,maxit);
itGMRESP = itGMRESP(2);
rrGMRESP = rrGMRESP(:)./rrGMRESP(1);
rrGMRESP = rrGMRESP(1:itout:end);
iterGMRESP = 0:itout:itout*size(rrGMRESP,1)-1;
rrGMRESP = [iterGMRESP' rrGMRESP];
xGMRESP = M\xGMRESP;
disp(['|             converges in ' num2str(itGMRESP) ' iterations']);



% Compute GMRES without prec
disp(['| GMRES...']);
[xGMRES, ~, ~, itGMRES, rrGMRES] = gmres(A,b,[],tol,maxit);
itGMRES = itGMRES(2);
rrGMRES = rrGMRES(:)./rrGMRES(1);
rrGMRES = rrGMRES(1:itout:end);
iterGMRES = 0:itout:itout*size(rrGMRES,1)-1;
rrGMRES = [iterGMRES' rrGMRES];
disp(['|             converges in ' num2str(itGMRES) ' iterations']);



%%%%%%%%%%% Deflation %%%%%%%%%%%



disp(['| Compute eigenvectors to deflate...']);
[eigvec,nbEigVec] = computeEigVec2D_cavity(mesh, dofm, nbEigVec,eigvecToDeflate);
% [eigvec,~] = eigs(A,nbEigVec,'smallestabs');
disp(['|             Done']);

[P,Q] = computeDefOp(nbEigVec, eigvec, A);


% disp(['| Compute eigenvalues of A*(P+Q)...']);
% [~, evdef] = eigs(A*(P+Q),50-nbEigVec,'smallestabs');
% evdef = diag(evdef);
% evdef = [real(evdef) imag(evdef)];
% % csvwrite(["output/evdef.csv"],evdef);
% disp(['|             Done']);



%%% No preconditioner :

% Compute GMRES with DEF1 and closest eigvec : P*A*x = P*b
% disp(['| Compute GMRES with DEF1 and closest eigvec...']);
% [xD, ~, ~, itD, rrD] = gmres(P*A,P*b,[],tol,maxit);
% itD = itD(2);
% rrD = rrD(:)./rrD(1);
% rrD = rrD(1:itout:end);
% iterD = 0:itout:itout*size(rrD,1)-1;
% rrD = [iterD' rrD];
% xD = P'*xD + Q*b;
% disp(['|             Done']);

% Compute GMRES with DEF1 and closest eigvec : A*P*u = b, x = P'*P*u + Q*b
% disp(['| GMRES with DEF1...']);
% [uD, ~, ~, itD, rrD] = gmres(A*P,b,[],tol,maxit);
% itD = itD(2);
% rrD = rrD(:)./rrD(1);
% rrD = rrD(1:itout:end);
% iterD = 0:itout:itout*size(rrD,1)-1;
% rrD = [iterD' rrD];
% xD = P'*P*uD + Q*b;
% disp(['|             converges in ' num2str(itD) ' iterations']);


% Compute GMRES with ADEF1 and closest eigvec : A*(P+Q)*u = b, x = (P+Q)*u
disp(['| GMRES with ADEF1...']);
[uAD, ~, ~, itAD, rrAD] = gmres(A*(P+Q),b,[],tol,maxit);
itAD = itAD(2);
rrAD = rrAD(:)./rrAD(1);
rrAD = rrAD(1:itout:end);
iterAD = 0:itout:itout*size(rrAD,1)-1;
rrAD = [iterAD' rrAD];
xAD = (P+Q)*uAD;
disp(['|             converges in ' num2str(itAD) ' iterations']);




%%% Add preconditioner :

MinvP = M\P;

% Compute GMRES with DEF1 and closest eigvec and prec : P*A/M*u = P*b, x = M\u
% disp(['| Compute GMRES with DEF1 and closest eigvec and prec...']);
% [xPD, ~, ~, itPD, rrPD] = gmres(P*A/M,P*b,[],tol,maxit);
% itPD = itPD(2);
% rrPD = rrPD(:)./rrPD(1);
% rrPD = rrPD(1:itout:end);
% iterPD = 0:itout:itout*size(rrPD,1)-1;
% rrPD = [iterPD' rrPD];
% xPD = P'/M*xPD + Q*b;
% disp(['|             Done']);


% Compute GMRES with DEF1 and closest eigvec and prec : A/M*P*u = b, x = P'\M*P*u + Q*b
% disp(['| Preconditoned GMRES with DEF1...']);
% [uPD, ~, ~, itPD, rrPD] = gmres(AMinv*P,b,[],tol,maxit);
% itPD = itPD(2);
% rrPD = rrPD(:)./rrPD(1);
% rrPD = rrPD(1:itout:end);
% iterPD = 0:itout:itout*size(rrPD,1)-1;
% rrPD = [iterPD' rrPD];
% xPD = P'*MinvP*uPD + Q*b;
% disp(['|             converges in ' num2str(itPD) ' iterations']);


% Compute GMRES with ADEF1 and closest eigvec and prec : A*(M\P+Q)*u = b, x = (M\P+Q)*u
disp(['| Preconditoned  GMRES with ADEF1...']);
[uPAD, ~, ~, itPAD, rrPAD] = gmres(A*(MinvP+Q),b,[],tol,maxit);
itPAD = itPAD(2);
rrPAD = rrPAD(:)./rrPAD(1);
rrPAD = rrPAD(1:itout:end);
iterPAD = 0:itout:itout*size(rrPAD,1)-1;
rrPAD = [iterPAD' rrPAD];
xPAD = (MinvP+Q)*uPAD;
disp(['|             converges in ' num2str(itPAD) ' iterations']);


%%% Save results %%%

folder = "output/freq_"+num2str(k);
if ~exist(folder, 'dir')
    mkdir(folder);
end

% for i=1:nbEigVec
%     writeField2D(dofm, mesh, eigvec(:,i), folder+"/eigvec"+num2str(i)+".pos", "eigvec"+num2str(i));
% end

writeField2D(dofm, mesh, xGMRESP, folder+"/solGMRESP.pos", "solGMRESP");
writeField2D(dofm, mesh, xGMRES, folder+"/solGMRES.pos", "solGMRES");
% writeField2D(dofm, mesh, xD, folder+"/solDef.pos", "solDef");
writeField2D(dofm, mesh, xAD, folder+"/solADef.pos", "solADef");
% writeField2D(dofm, mesh, xPD, folder+"/solDefP.pos", "solDefP");
writeField2D(dofm, mesh, xPAD, folder+"/solADefP.pos", "solADefP");

csvwrite([folder+"/rrGMRES.csv"],rrGMRES);
csvwrite([folder+"/rrGMRESP.csv"],rrGMRESP);
% csvwrite([folder+"/rrD.csv"],rrD);
csvwrite([folder+"/rrAD.csv"],rrAD);
% csvwrite([folder+"/rrPD.csv"],rrPD);
csvwrite([folder+"/rrPAD.csv"],rrPAD);

% return;

% csvwrite([folder+"/evA.csv"],evA);
% csvwrite([folder+"/evdef.csv"],evdef);

% it = [ itGMRES itGMRESP itD itAD itPD itPAD];
it = [ itGMRES itGMRESP itAD itPAD];

% csvwrite([folder+"/it.csv"],it');

%%% Plot results %%%

green = [0.4660 0.6740 0.1880];
magenta = [0.4940 0.1840 0.5560];
orange = [0.9290 0.6940 0.1250];
cyan = [0.3010 0.7450 0.9330];

%%% Spectra


% figure;
% hold on;
% s1 = scatter(real(evA),imag(evA),100,'DisplayName','Eigenvalues of A');
% s1.Marker = '+';
% s1.MarkerEdgeColor = 'b';
% 
% s2 = scatter(real(evAMinv),imag(evAMinv),100, 'DisplayName','Eigenvalues of P\A');
% s2.Marker = 'x';
% s2.MarkerEdgeColor = 'k';
% 
% s3 = scatter(real(evdef),imag(evdef),100, 'DisplayName','Eigenvalues of (P+Q)*A');
% s3.Marker = 'o';
% s3.MarkerEdgeColor = 'r';
% 
% legend('Location', 'southwest', 'fontsize', 15)
% 
% grid on; box on;
% title(['Spectra of $A$ - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)], 'interpreter', 'latex', 'fontsize', 20)
% % xlim([-0.04 0.06]);
% % hold on; plot(fovals(A));
% set(0,'DefaultFigureWindowStyle','docked')



%%% Residuals

maxIt = max(it);
minIt = min(it);

disp(['GMRES: ' num2str(itGMRES)]);
disp(['GMRES and Shift: ' num2str(itGMRESP)]);
% disp(['DEF1 : ' num2str(itD)]);
disp(['ADEF1 : ' num2str(itAD)]);
% disp(['DEF1 and Shift : ' num2str(itPD)]);
disp(['ADEF1 and Shift : ' num2str(itPAD)]);


figure
hold on
set(0,'DefaultFigureWindowStyle','docked')

p1 = semilogy(rrGMRES(:,1),rrGMRES(:,2),'b-o','DisplayName','Relative residual','linewidth', 2,'markersize', 10);
p2 = semilogy(rrGMRESP(:,1),rrGMRESP(:,2),'r-o','DisplayName','Relative residual with shift','linewidth', 2,'markersize', 10);
% p3 = semilogy(rrD(:,1),rrD(:,2),'g-o','DisplayName','Relative residual with DEF1','linewidth', 2,'markersize', 10);
% p3.Color = green;
p4 = semilogy(rrAD(:,1),rrAD(:,2),'c-o','DisplayName','Relative residual with ADEF1','linewidth', 2,'markersize', 10);
p4.Color = cyan;
% p5 = semilogy(rrPD(:,1),rrPD(:,2),'m-o','DisplayName','Relative residual with DEF1 and shift','linewidth', 2,'markersize', 10);
% p5.Color = magenta;
p6  = semilogy(rrPAD(:,1),rrPAD(:,2),'y-o','DisplayName','Relative residual with ADEF1 and shift','linewidth', 2,'markersize', 10);
p6.Color = orange;
% plot([0 maxit],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (direct)');

set(gca, 'YScale', 'log')
box on
grid on
xlim([0 maxIt+1]);
ylim auto;
title(['CG - ' benchmark ' - GMRES - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree) ' - nbEigvec=' num2str(nbEigVec)], 'interpreter', 'latex', 'fontsize', 20)
xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
legend('Location', 'southwest', 'fontsize', 15)