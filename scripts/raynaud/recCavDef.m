
clear all;
%close all;


disp(['---------------------------------------------------------']);
disp(['WARNING : deflation for Neumann BC not implemented']);
disp(['---------------------------------------------------------']);


global k h

N=15;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);

computeSolNum2D = @computeSolNum2D_CG;

% Setup benchmark and parameters
benchmark = 'scattering_rec';
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
        k = 23.597;
%         k = 10;
        h = 1/20;
        tol = 1e-6; maxit = 5000; itout = 10;
        LdomX = 0.95;
        LdomY = 0.5;
        LpmlX = 0.2;
        LpmlY = 0.2;
    case 'waveguide'
        k = 6*pi;
        h = 1/8;
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 3; % P1
PREC = 1; % for preconditioner
nbEigVec=1;

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

disp(['| Compute system and deflation subspace...']);

[~, sysA] = computeSolNum2D(mesh, dofm, PREC);

A = sysA.matA;
M = sysA.matP;
b = sysA.rhsA;
AMinv = A/M;

[eigvec,~] = eigs(A,nbEigVec,'smallestabs');

[P,Q] = computeDefOp(nbEigVec, eigvec, A);

MinvP = M\P;

if maxit > size(A,2)
    maxit = size(A,2);
end


disp(['|             Done']);

% disp(['| Compute eigenvalues of A...']);
% 
% [~, evA] = eigs(A,size(A,2),'smallestabs');
% evA = diag(evA);
% evA = [real(evA) imag(evA)];
% csvwrite(["output/evA.csv"],evA);

% disp(['|             Done']);


% disp(['| Compute eigenvalues of A/M...']);
% [~, evAMinv] = eigs(AMinv,50,'smallestabs');
% evAMinv = diag(evAMinv);
% evAMinv = [real(evAMinv) imag(evAMinv)];
% csvwrite(["output/spectrum_prec.csv"],evAMinv);
% disp(['|             Done']);


% disp(['| Compute eigenvalues of A*(P+Q)...']);
% [~, evD] = eigs(A*(P+Q),50-nbEigVec,'smallestabs');
% evD = diag(evD);
% evD = [real(evD) imag(evD)];
% disp(['|             Done']);

% disp(['| Compute eigenvalues of A*(M\P+Q)...']);
% [~, evPD] = eigs(A*(MinvP+Q),50-nbEigVec,'smallestabs');
% evPD = diag(evPD);
% evPD = [real(evPD) imag(evPD)];
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


%%% No preconditioner :


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

% compute incident field
disp(['Compute incident field']);

solInc = -computeSolProjL2_2D_CG(mesh, dofm);

%%% Save results %%%

folder = "output/freq_"+num2str(k);
if ~exist(folder, 'dir')
    mkdir(folder);
end

% for i=1:nbEigVec
%     writeField2D(dofm, mesh, eigvec(:,i), folder+"/eigvec"+num2str(i)+".pos", "eigvec"+num2str(i));
% end

global WRITE_FIELD_ABSOLUTE
WRITE_FIELD_ABSOLUTE = 1;

writeField2D(dofm, mesh, xGMRESP+solInc, folder+"/solGMRESP.pos", "solGMRESP");
writeField2D(dofm, mesh, xGMRES+solInc, folder+"/solGMRES.pos", "solGMRES");
writeField2D(dofm, mesh, xAD+solInc, folder+"/solADef.pos", "solADef");
writeField2D(dofm, mesh, xPAD+solInc, folder+"/solADefP.pos", "solADefP");

csvwrite([folder+"/rrGMRES.csv"],rrGMRES);
csvwrite([folder+"/rrGMRESP.csv"],rrGMRESP);
csvwrite([folder+"/rrAD.csv"],rrAD);
csvwrite([folder+"/rrPAD.csv"],rrPAD);


% csvwrite([folder+"/evA.csv"],evA);
% csvwrite([folder+"/evAMinv.csv"],evAMinv);
% csvwrite([folder+"/evD.csv"],evD);
% csvwrite([folder+"/evPD.csv"],evPD);

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
% s3 = scatter(real(evD),imag(evD),100, 'DisplayName','Eigenvalues of (P+Q)*A');
% s3.Marker = 'o';
% s3.MarkerEdgeColor = 'r';
% 
% s4 = scatter(real(evPD),imag(evPD),100, 'DisplayName','Eigenvalues of (M\P+Q)*A');
% s4.Marker = 'o';
% s4.MarkerEdgeColor = 'g';
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
disp(['ADEF1 : ' num2str(itAD)]);
disp(['ADEF1 and Shift : ' num2str(itPAD)]);


figure
hold on
set(0,'DefaultFigureWindowStyle','docked')

p1 = semilogy(rrGMRES(:,1),rrGMRES(:,2),'b-o','DisplayName','Relative residual','linewidth', 2,'markersize', 10);
p2 = semilogy(rrGMRESP(:,1),rrGMRESP(:,2),'r-o','DisplayName','Relative residual with shift','linewidth', 2,'markersize', 10);
p3 = semilogy(rrAD(:,1),rrAD(:,2),'c-o','DisplayName','Relative residual with ADEF1','linewidth', 2,'markersize', 10);
p3.Color = cyan;
p4  = semilogy(rrPAD(:,1),rrPAD(:,2),'y-o','DisplayName','Relative residual with ADEF1 and shift','linewidth', 2,'markersize', 10);
p4.Color = orange;
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