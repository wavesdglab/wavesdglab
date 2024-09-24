
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
h = 1/64;
tol = 1e-8; maxit = 2000; itout =4;
L = 1;

degree = 1;
PREC = 1;

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
prec = sysA.matP;
b = sysA.rhsA;

[~, eigenvalA] = eigs(A,30,'smallestabs');
eigenvalA = diag(eigenvalA);

%%%%%%%%%%% No deflation %%%%%%%%%%%

% Compute GMRES with prec
[resGMRESP, errGMRESP, itGMRESP, ~, ~, xGMRESP] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);


% Compute GMRES without prec
sysA.matP = speye(size(A,1));
[resGMRES, errGMRES, itGMRES, ~, ~, xGMRES] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

%%%%%%%%%%% Closest deflated eigvec %%%%%%%%%%%

nbEigVec=22;
[eigenvecC,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,"closestEigvec");
[eigenvecC,~] = eigs(A,nbEigVec,'smallestabs');


[P,Q] = computeDefOp(nbEigVec, eigenvecC, A);

% [~, eigenvalPA_C] = eigs(P*A,30,'smallestabs');
% eigenvalPA_C = diag(eigenvalPA_C);

[~, eigenvalPAQ_C] = eigs(P*A+Q,30,'smallestabs');
eigenvalPAQ_C = diag(eigenvalPAQ_C);

%%% No preconditioner :

% Compute GMRES with DEF1 and closest eigvec
sysA.matA = P*A;
sysA.rhsA = P*b;
[resDef1C, errDef1C, itDef1C, ~, ~, xDef1C] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES with ADEF1 and closest eigvec
sysA.matA = (P+Q)*A;
sysA.rhsA = (P+Q)*b;
[resADef1C, errADef1C, itADef1C, ~, ~, xADef1C] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

%%% Add preconditioner :

% Compute GMRES with DEF1 and closest eigvec and prec
sysA.matP = prec;
sysA.matA = P*A;
sysA.rhsA = P*b;
[resDef1PC, errDef1PC, itDef1PC, ~, ~, xDef1PC] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES with ADEF1M and closest eigvec and prec
sysA.matA = (P+Q)*A;
sysA.rhsA = (P+Q)*b;
[resADef1PC, errADef1PC, itADef1PC, ~, ~, xADef1PC] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);


% Compute GMRES with ADEF1 and closest eigvec and prec
sysA.matP = speye(size(A,1));
PP = prec\P;
sysA.matA = (PP+Q)*A;
sysA.rhsA = (PP+Q)*b;
[resADef1MPC, errADef1MPC, itADef1MPC, ~, ~, xADef1MPC] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);


%%%%%%%%%%% First deflated eigvec %%%%%%%%%%%

[eigenvecF,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,"firstEigvec");
[eigenvecF,~] = eigs(A,nbEigVec,'smallestreal');

[P,Q] = computeDefOp(nbEigVec, eigenvecF, A);

% [~, eigenvalPA_F] = eigs(P*A,30,'smallestabs');
% eigenvalPA_F = diag(eigenvalPA_F);

[~, eigenvalPAQ_F] = eigs(P*A+Q,30,'smallestabs');
eigenvalPAQ_F = diag(eigenvalPAQ_F);

%%% No preconditioner :

% Compute GMRES with DEF1 and first eigvec
sysA.matA = P*A;
sysA.rhsA = P*b;
[resDef1F, errDef1F, itDef1F, ~, ~, xDef1F] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES with ADEF1 and first eigvec
sysA.matA = (P+Q)*A;
sysA.rhsA = (P+Q)*b;
[resADef1F, errADef1F, itADef1F, ~, ~, xADef1F] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

%%% Add preconditioner :

% Compute GMRES with DEF1 and first eigvec and prec
sysA.matP = prec;
sysA.matA = P*A;
sysA.rhsA = P*b;
[resDef1PF, errDef1PF, itDef1PF, ~, ~, xDef1PF] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES with ADEF1M and first eigvec and prec
sysA.matA = (P+Q)*A;
sysA.rhsA = (P+Q)*b;
[resADef1PF, errADef1PF, itADef1PF, ~, ~, xADef1PF] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

% Compute GMRES with ADEF1 and first eigvec and prec
sysA.matP = speye(size(A,1));
PP = prec\P;
sysA.matA = (PP+Q)*A;
sysA.rhsA = (PP+Q)*b;
[resADef1MPF, errADef1MPF, itADef1MPF, ~, ~, xADef1MPF] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);


%%% Save results %%%
% 
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
% 
% folder = "output/eigvec"+num2str(nbEigVec);
% if ~exist(folder, 'dir')
%     mkdir(folder);
% end
% csvwrite([folder+"/resGMRES.csv"],resGMRES);
% csvwrite([folder+"/resGMRESP.csv"],resGMRESP);
% csvwrite([folder+"/resDef1C.csv"],resDef1C);
% csvwrite([folder+"/resADef1C.csv"],resADef1C);
% csvwrite([folder+"/resDef1PC.csv"],resDef1PC);
% csvwrite([folder+"/resADef1PC.csv"],resADef1PC);
% csvwrite([folder+"/resADef1MPC.csv"],resADef1MPC);
% csvwrite([folder+"/resDef1F.csv"],resDef1F);
% csvwrite([folder+"/resADef1F.csv"],resADef1F);
% csvwrite([folder+"/resDef1PF.csv"],resDef1PF);
% csvwrite([folder+"/resADef1PF.csv"],resADef1PF);
% csvwrite([folder+"/resADef1MPF.csv"],resADef1MPF);
%
% csvwrite([folder+"/errGMRES.csv"],errGMRES);
% csvwrite([folder+"/errGMRESP.csv"],errGMRESP);
% csvwrite([folder+"/errDef1C.csv"],errDef1C);
% csvwrite([folder+"/errADef1C.csv"],errADef1C);
% csvwrite([folder+"/errDef1PC.csv"],errDef1PC);
% csvwrite([folder+"/errADef1PC.csv"],errADef1PC);
% csvwrite([folder+"/errADef1MPC.csv"],errADef1MPC);
% csvwrite([folder+"/errDef1F.csv"],errDef1F);
% csvwrite([folder+"/errADef1F.csv"],errADef1F);
% csvwrite([folder+"/errDef1PF.csv"],errDef1PF);
% csvwrite([folder+"/errADef1PF.csv"],errADef1PF);
% csvwrite([folder+"/errADef1MPF.csv"],errADef1MPF);
% 
% csvwrite([folder+"/eigenvalA.csv"],eigenvalA);
% % csvwrite([folder+"/eigenvalPA_C.csv"],eigenvalPA_C);
% csvwrite([folder+"/eigenvalPAQ_C.csv"],eigenvalPAQ_C);
% % csvwrite([folder+"/eigenvalPA_F.csv"],eigenvalPA_F);
% csvwrite([folder+"/eigenvalPAQ_F.csv"],eigenvalPAQ_F);

%%% Plot results %%%

green = [0.4660 0.6740 0.1880];
magenta = [0.4940 0.1840 0.5560];
orange = [0.9290 0.6940 0.1250];
cyan = [0.3010 0.7450 0.9330];

%%% Spectra


figure;
hold on;
s1 = scatter(real(eigenvalA),imag(eigenvalA),100,'DisplayName','Eigenvalues of A');
s1.Marker = '+';
s1.MarkerEdgeColor = 'b';

% s2 = scatter(real(eigenvalPA_C),imag(eigenvalPA_C),100, 'DisplayName','Eigenvalues of PA, closest eigvec');
% s2.Marker = 'o';
% s2.MarkerEdgeColor = green;

s3 = scatter(real(eigenvalPAQ_C),imag(eigenvalPAQ_C),100, 'DisplayName','Eigenvalues of PA+Q, closest eigvec');
s3.Marker = 'o';
s3.MarkerEdgeColor = green;

% s4 = scatter(real(eigenvalPA_F),imag(eigenvalPA_F),100, 'DisplayName','Eigenvalues of PA, first eigvec');
% s4.Marker = 'x';
% s4.MarkerEdgeColor = green;

s5 = scatter(real(eigenvalPAQ_F),imag(eigenvalPAQ_F),100, 'DisplayName','Eigenvalues of PA+Q, first eigvec');
s5.Marker = 'x';
s5.MarkerEdgeColor = green;

legend('Location', 'southwest', 'fontsize', 15)

grid on; box on;
title(['Spectra of $A$ and $PA+Q$ - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree) ' - nbEigvec=' num2str(nbEigVec)], 'interpreter', 'latex', 'fontsize', 20)
xlim([-0.04 0.06]);
set(0,'DefaultFigureWindowStyle','docked')



%%% Residuals

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

p1 = semilogy(0:itout:maxit,resGMRES,'b-o','DisplayName','Relative residual','linewidth', 2,'markersize', 10);
p2 = semilogy(0:itout:maxit,resGMRESP,'r-o','DisplayName','Relative residual with shift','linewidth', 2,'markersize', 10);
p3 = semilogy(0:itout:maxit,resDef1C,'g-o','DisplayName','Relative residual with DEF1 closest','linewidth', 2,'markersize', 10);
p3.Color = green;
p4 = semilogy(0:itout:maxit,resADef1C,'c-o','DisplayName','Relative residual with ADEF1 closest','linewidth', 2,'markersize', 10);
p4.Color = cyan;
p5 = semilogy(0:itout:maxit,resDef1PC,'m-o','DisplayName','Relative residual with DEF1 closest and shift','linewidth', 2,'markersize', 10);
p5.Color = magenta;
p6  = semilogy(0:itout:maxit,resADef1PC,'y-o','DisplayName','Relative residual with ADEF1 closest and shift','linewidth', 2,'markersize', 10);
p6.Color = orange;
p7 = semilogy(0:itout:maxit,resADef1MPC,'k-o','DisplayName','Relative residual with ADEF1M closest and shift','linewidth', 2,'markersize', 10);
p8 = semilogy(0:itout:maxit,resDef1F,'g-x','DisplayName','Relative residual with DEF1 first','linewidth', 2,'markersize', 10);
p8.Color = green;
p9 = semilogy(0:itout:maxit,resADef1F,'c-x','DisplayName','Relative residual with ADEF1 first','linewidth', 2,'markersize', 10);
p9.Color = cyan;
p10 = semilogy(0:itout:maxit,resDef1PF,'m-x','DisplayName','Relative residual with DEF1 first and shift','linewidth', 2,'markersize', 10);
p10.Color = magenta;
p11 = semilogy(0:itout:maxit,resADef1PF,'y-x','DisplayName','Relative residual with ADEF1 first and shift','linewidth', 2,'markersize', 10);
p11.Color = orange;
p12 = semilogy(0:itout:maxit,resADef1MPF,'k-x','DisplayName','Relative residual with ADEF1M first and shift','linewidth', 2,'markersize', 10);

p13 = semilogy(0:itout:maxit,errGMRES,'b--','DisplayName','Relative error','linewidth', 2,'markersize', 10);
p14 = semilogy(0:itout:maxit,errGMRESP,'r--','DisplayName','Relative error with shift','linewidth', 2,'markersize', 10);
p15 = semilogy(0:itout:maxit,errDef1C,'g--','DisplayName','Relative error with DEF1 closest','linewidth', 2,'markersize', 10);
p15.Color = green;
p16 = semilogy(0:itout:maxit,errADef1C,'c--','DisplayName','Relative error with ADEF1 closest','linewidth', 2,'markersize', 10);
p16.Color = cyan;
p17 = semilogy(0:itout:maxit,errDef1PC,'m--','DisplayName','Relative error with DEF1 closest and shift','linewidth', 2,'markersize', 10);
p17.Color = magenta;
p18 = semilogy(0:itout:maxit,errADef1PC,'y--','DisplayName','Relative error with ADEF1 closest and shift','linewidth', 2,'markersize', 10);
p18.Color = orange;
p19 = semilogy(0:itout:maxit,errADef1MPC,'k--','DisplayName','Relative error with ADEF1M closest and shift','linewidth', 2,'markersize', 10);
p20 = semilogy(0:itout:maxit,errDef1F,'g--','DisplayName','Relative error with DEF1 first','linewidth', 2,'markersize', 10);
p20.Color = green;
p21 = semilogy(0:itout:maxit,errADef1F,'c--','DisplayName','Relative error with ADEF1 first','linewidth', 2,'markersize', 10);
p21.Color = cyan;
p22 = semilogy(0:itout:maxit,errDef1PF,'m--','DisplayName','Relative error with DEF1 first and shift','linewidth', 2,'markersize', 10);
p22.Color = magenta;
p23 = semilogy(0:itout:maxit,errADef1PF,'y--','DisplayName','Relative error with ADEF1 first and shift','linewidth', 2,'markersize', 10);
p23.Color = orange;
p24 = semilogy(0:itout:maxit,errADef1MPF,'k--','DisplayName','Relative error with ADEF1M first and shift','linewidth', 2,'markersize', 10);



set(gca, 'YScale', 'log')
box on
grid on
xlim([0 maxIt+1]);
ylim auto;
title(['CG - ' benchmark ' - GMRES - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree) ' - nbEigvec=' num2str(nbEigVec)], 'interpreter', 'latex', 'fontsize', 20)
xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
legend('Location', 'southwest', 'fontsize', 15)