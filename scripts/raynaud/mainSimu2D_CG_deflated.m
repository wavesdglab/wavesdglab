
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
        k = 2.001*sqrt(2)*pi;
%         k = 5.877;
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
degree = 1; % P1
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
prec = sysA.matP;
b = sysA.rhsA;

[~, eigenvalA] = eigs(A,10,'smallestabs');
eigenvalA = diag(eigenvalA);
figure;
hold off; scatter(real(eigenvalA),imag(eigenvalA));
grid on; box on;
title(['Spectum of A'], 'interpreter', 'latex', 'fontsize', 20)
xlim([-0.02 0.02]);
set(0,'DefaultFigureWindowStyle','docked')

%%%%%%%%%%% No deflation %%%%%%%%%%%

% Compute GMRES with prec
[resVecGMRESP, ~, itGMRESP, ~, ~, xGMRESP] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES without prec
sysA.matP = speye(size(A,1));
[resVecGMRES, ~, itGMRES, ~, ~, xGMRES] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

%%%%%%%%%%% Closest deflated eigvec %%%%%%%%%%%

nbEigVec=5;
[eigenvecC,nbEigVec] = computeEigVec2D_cavity(mesh, dofm, nbEigVec,"closesteigvec");

[P,Q] = computeDefOp(nbEigVec, eigenvecC, A);

[~, eigenvalPA] = eigs(P*A,10,'smallestabs');
eigenvalPA = diag(eigenvalPA);
figure;
hold off; scatter(real(eigenvalPA),imag(eigenvalPA));
grid on; box on;
title(['Spectum of PA, closest deflated eigvec'], 'interpreter', 'latex', 'fontsize', 20)
xlim([-0.02 0.02]);
set(0,'DefaultFigureWindowStyle','docked')

[~, eigenvalPAQ] = eigs(P*A+Q,10,'smallestabs');
eigenvalPAQ = diag(eigenvalPAQ);
figure;
hold off; scatter(real(eigenvalPAQ),imag(eigenvalPAQ));
grid on; box on;
title(['Spectum of PA+Q, closest deflated eigvec'], 'interpreter', 'latex', 'fontsize', 20)
xlim([-0.02 0.02]);
set(0,'DefaultFigureWindowStyle','docked')

%%% No preconditioner :

% Compute GMRES with DEF1 and closest eigvec
sysA.matA = P*A;
sysA.rhsA = P*b;
[resVecDef1C, ~, itDef1C, ~, ~, xDef1C] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES with ADEF1 and closest eigvec
sysA.matA = (P+Q)*A;
sysA.rhsA = (P+Q)*b;
[resVecADef1C, ~, itADef1C, ~, ~, xADef1C] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

%%% Add preconditioner :

% Compute GMRES with DEF1 and closest eigvec and prec
sysA.matP = prec;
sysA.matA = P*A;
sysA.rhsA = P*b;
[resVecDef1PC, ~, itDef1PC, ~, ~, xDef1PC] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES with ADEF1M and closest eigvec and prec
sysA.matA = (P+Q)*A;
sysA.rhsA = (P+Q)*b;
[resVecADef1PC, ~, itADef1PC, ~, ~, xADef1PC] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES with ADEF1 and closest eigvec and prec
sysA.matP = speye(size(A,1));
PP = prec\P;
sysA.matA = (PP+Q)*A;
sysA.rhsA = (PP+Q)*b;
[resVecADef1MPC, ~, itADef1MPC, ~, ~, xADef1MPC] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);


%%%%%%%%%%% First deflated eigvec %%%%%%%%%%%

[eigenvecF,nbEigVec] = computeEigVec2D_cavity(mesh, dofm, nbEigVec,"firsteigvec");

[P,Q] = computeDefOp(nbEigVec, eigenvecF, A);

[~, eigenvalPA] = eigs(P*A,10,'smallestabs');
eigenvalPA = diag(eigenvalPA);
figure;
hold off; scatter(real(eigenvalPA),imag(eigenvalPA));
grid on; box on;
title(['Spectum of PA, first deflated eigvec'], 'interpreter', 'latex', 'fontsize', 20)
xlim([-0.02 0.02]);
set(0,'DefaultFigureWindowStyle','docked')

[~, eigenvalPAQ] = eigs(P*A+Q,10,'smallestabs');
eigenvalPAQ = diag(eigenvalPAQ);
figure;
hold off; scatter(real(eigenvalPAQ),imag(eigenvalPAQ));
grid on; box on;
title(['Spectum of PA+Q, first deflated eigvec'], 'interpreter', 'latex', 'fontsize', 20)
xlim([-0.02 0.02]);
set(0,'DefaultFigureWindowStyle','docked')


%%% No preconditioner :

% Compute GMRES with DEF1 and first eigvec
sysA.matA = P*A;
sysA.rhsA = P*b;
[resVecDef1F, ~, itDef1F, ~, ~, xDef1F] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES with ADEF1 and first eigvec
sysA.matA = (P+Q)*A;
sysA.rhsA = (P+Q)*b;
[resVecADef1F, ~, itADef1F, ~, ~, xADef1F] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

%%% Add preconditioner :

% Compute GMRES with DEF1 and first eigvec and prec
sysA.matP = prec;
sysA.matA = P*A;
sysA.rhsA = P*b;
[resVecDef1PF, ~, itDef1PF, ~, ~, xDef1PF] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES with ADEF1M and first eigvec and prec
sysA.matA = (P+Q)*A;
sysA.rhsA = (P+Q)*b;
[resVecADef1PF, ~, itADef1PF, ~, ~, xADef1PF] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES with ADEF1 and first eigvec and prec
sysA.matP = speye(size(A,1));
PP = prec\P;
sysA.matA = (PP+Q)*A;
sysA.rhsA = (PP+Q)*b;
[resVecADef1MPF, ~, itADef1MPF, ~, ~, xADef1MPF] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);


%%% Save results

% writeField2D(dofm, mesh, xGMRESP, 'output/solGMRESP.pos', "solGMRESP");
% writeField2D(dofm, mesh, xGMRES, 'output/solGMRES.pos', "solGMRES");
% writeField2D(dofm, mesh, xDef1C, 'output/solDef1C.pos', "solDef1C");
% writeField2D(dofm, mesh, xADef1C, 'output/solADef1C.pos', "solADef1C");
% writeField2D(dofm, mesh, xDef1PC, 'output/solDef1PC.pos', "solDef1PC");
% writeField2D(dofm, mesh, xADef1PC, 'output/solADef1PC.pos', "solADef1PC");
% writeField2D(dofm, mesh, xADef1MPC, 'output/solADef1MPC.pos', "solADef1MPC");
% writeField2D(dofm, mesh, xDef1F, 'output/solDef1F.pos', "solDef1F");
% writeField2D(dofm, mesh, xADef1F, 'output/solADef1F.pos', "solADef1F");
% writeField2D(dofm, mesh, xDef1PF, 'output/solDef1PF.pos', "solDef1PF");
% writeField2D(dofm, mesh, xADef1PF, 'output/solADef1PF.pos', "solADef1PF");
% writeField2D(dofm, mesh, xADef1MPF, 'output/solADef1MPF.pos', "solADef1MPF");
% system('gmsh output/solGMRESP.pos output/solGMRES.pos output/solDef1C.pos output/solADef1C.pos output/solDef1PC.pos output/solADef1PC.pos output/solADef1MPC.pos output/solDef1F.pos output/solADef1F.pos output/solDef1PF.pos output/solADef1PF.pos output/solADef1MPF.pos&');

folder = "output/eigvec"+num2str(nbEigVec);
csvwrite([folder+"/resVecGMRES.csv"],resVecGMRES);
csvwrite([folder+"/resVecGMRESP.csv"],resVecGMRESP);
csvwrite([folder+"/resVecDef1C.csv"],resVecDef1C);
csvwrite([folder+"/resVecADef1C.csv"],resVecADef1C);
csvwrite([folder+"/resVecDef1PC.csv"],resVecDef1PC);
csvwrite([folder+"/resVecADef1PC.csv"],resVecADef1PC);
csvwrite([folder+"/resVecADef1MPC.csv"],resVecADef1MPC);
csvwrite([folder+"/resVecDef1F.csv"],resVecDef1F);
csvwrite([folder+"/resVecADef1F.csv"],resVecADef1F);
csvwrite([folder+"/resVecDef1PF.csv"],resVecDef1PF);
csvwrite([folder+"/resVecADef1PF.csv"],resVecADef1PF);
csvwrite([folder+"/resVecADef1MPF.csv"],resVecADef1MPF);


%%% Plot results

it = [itGMRESP itGMRES itDef1C itADef1C itDef1PC itADef1PC itADef1MPC itDef1F itADef1F itDef1PF itADef1PF itADef1MPF];

maxIt = max(it);
minIt = min(it);

disp(['GMRES: ' num2str(itGMRES)]);
disp(['GMRES and Shift: ' num2str(itGMRESP)]);
disp(['DEF1 Closest: ' num2str(itDef1C)]);
disp(['DEF1 First: ' num2str(itDef1F)]);
disp(['ADEF1 Closest: ' num2str(itADef1C)]);
disp(['ADEF1 First: ' num2str(itADef1F)]);
disp(['DEF1 and Shift Closest: ' num2str(itDef1PC)]);
disp(['DEF1 and Shift First: ' num2str(itDef1PF)]);
disp(['ADEF1 and Shift Closest: ' num2str(itADef1PC)]);
disp(['ADEF1 and Shift First: ' num2str(itADef1PF)]);
disp(['ADEF1M and Shift Closest: ' num2str(itADef1MPC)]);
disp(['ADEF1M and Shift First: ' num2str(itADef1MPF)]);

disp(['Difference: ' num2str(100*(maxIt-minIt)/maxIt) '%']);


figure
hold on
set(0,'DefaultFigureWindowStyle','docked')

green = [0.4660 0.6740 0.1880];
magenta = [0.4940 0.1840 0.5560];
orange = [0.9290 0.6940 0.1250];
cyan = [0.3010 0.7450 0.9330];

p1 = semilogy(0:itout:maxit,resVecGMRES,'b-o','DisplayName','Relative residual','linewidth', 2,'markersize', 10);
p2 = semilogy(0:itout:maxit,resVecGMRESP,'r-o','DisplayName','Relative residual with shift','linewidth', 2,'markersize', 10);
p3 = semilogy(0:itout:maxit,resVecDef1C,'g-o','DisplayName','Relative residual with DEF1 closest','linewidth', 2,'markersize', 10);
p3.Color = green;
p4 = semilogy(0:itout:maxit,resVecADef1C,'c-o','DisplayName','Relative residual with ADEF1 closest','linewidth', 2,'markersize', 10);
p4.Color = cyan;
p5 = semilogy(0:itout:maxit,resVecDef1PC,'m-o','DisplayName','Relative residual with DEF1 closest and shift','linewidth', 2,'markersize', 10);
p5.Color = magenta;
p6  = semilogy(0:itout:maxit,resVecADef1PC,'y-o','DisplayName','Relative residual with ADEF1 closest and shift','linewidth', 2,'markersize', 10);
p6.Color = orange;
p7 = semilogy(0:itout:maxit,resVecADef1MPC,'k-o','DisplayName','Relative residual with ADEF1M closest and shift','linewidth', 2,'markersize', 10);
p8 = semilogy(0:itout:maxit,resVecDef1F,'g-x','DisplayName','Relative residual with DEF1 first','linewidth', 2,'markersize', 10);
p8.Color = green;
p9 = semilogy(0:itout:maxit,resVecADef1F,'c-x','DisplayName','Relative residual with ADEF1 first','linewidth', 2,'markersize', 10);
p9.Color = cyan;
p10 = semilogy(0:itout:maxit,resVecDef1PF,'m-x','DisplayName','Relative residual with DEF1 first and shift','linewidth', 2,'markersize', 10);
p10.Color = magenta;
p11 = semilogy(0:itout:maxit,resVecADef1PF,'y-x','DisplayName','Relative residual with ADEF1 first and shift','linewidth', 2,'markersize', 10);
p11.Color = orange;
p12 = semilogy(0:itout:maxit,resVecADef1MPF,'k-x','DisplayName','Relative residual with ADEF1M first and shift','linewidth', 2,'markersize', 10);


set(gca, 'YScale', 'log')
box on
grid on
xlim([0 maxIt+1]);
ylim auto;
title(['CG - ' benchmark ' - GMRES - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)'], 'interpreter', 'latex', 'fontsize', 20)
xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
legend('Location', 'southwest', 'fontsize', 15)