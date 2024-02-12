close all;
clear all;

global h
global omega eta1 eta2 k1 k2 c1 c2 rho1 rho2
global rho c eta k

degree = 3;
BASIS = 1;
PREC = 1;
A = 1;              % order of numerical fluxes
B = 2;              % order of transmission variables
tol = 1e-100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLANE WAVE IN HOMOGENEOUS MEDIUM WITH EXACT DIRICHLET B.C.
iMax = 1000; iOut = 50;
benchmark = 'cavity_heterogeneous';

degree = 3; h = 1/10; omega = 15*pi;
rho1 = 1; c1 = 1; rho2 = 1; c2 = 1;
eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
run(benchmark,degree,BASIS,PREC,A,B,tol,iMax,iOut);

degree = 3; h = 1/10; omega = 15*pi;
rho1 = 1; c1 = 0.5; rho2 = 1; c2 = 0.5;
eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
run(benchmark,degree,BASIS,PREC,A,B,tol,iMax,iOut);

% PLANE WAVE IN INHOMOGENEOUS MEDIUM WITH EXACT DIRICHLET B.C.
iMax = 1000; iOut = 50;
benchmark = 'cavity_heterogeneous';

degree = 3; h = 1/10; omega = 15*pi;
rho1 = 1; c1 = 1; rho2 = 1; c2 = 0.8;
eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
run(benchmark,degree,BASIS,PREC,A,B,tol,iMax,iOut);

degree = 3; h = 1/10; omega = 15*pi;
rho1 = 1; c1 = 2; rho2 = 1; c2 = 1.6;
eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
run(benchmark,degree,BASIS,PREC,A,B,tol,iMax,iOut);

% BENCH WAVEGUIDE
% iMax = 4000; iOut = 200;
% benchmark = 'waveguide'; k = 6*pi; h = 1/8;
% run(benchmark,degree,tau,BASIS,PREC,tol,iMax,iOut);
% benchmark = 'waveguide'; k = 12*pi; h = 1/17;
% run(benchmark,degree,tau,BASIS,PREC,tol,iMax,iOut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,BASIS,PREC,A,B,tol,iMax,iOut)
global h
global omega eta1 eta2 k1 k2 c1 c2 rho1 rho2
global rho c eta k

% Build mesh and dofManager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);
setParameters(mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method CHDG (' benchmark ')']);
disp(['---------------------------------------------------------']);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_CHDG_ALL(mesh, dofm, PREC, A, B);
normErr= computeNormError2D_DG_ALL(mesh, dofm, solA);

solP = computeSolProjL2_2D_DG(mesh, dofm);
normProjErr = computeNormError2D_DG_ALL(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(normErr,'%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(normProjErr,'%1.2e')]);
disp(['---------------------------------------------------------']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver Richardson']);
alpha = 1;
[resRedVec, resPhyVec, error0] = solverRichardsonRedu_DG(mesh, dofm, sysA, tol, iMax, iOut, alpha, @computeNormError2D_DG);

iterVec = (0:iOut:iMax)';
errorRef = normErr*ones(size(error0));

rezu1 = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2 = [iterVec resRedVec, resPhyVec, error0, errorRef];
name = sprintf('output/historyRich_CHDG_%s_p%i_h%i_alpha%g.csv', benchmark, degree, h, alpha);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver CGNR']);
[resRedVec, resPhyVec, error1] = solverCGNRredu_DG(mesh, dofm, sysA, tol, iMax, iOut, @computeNormError2D_DG);

iterVec = (0:iOut:iMax)';
errorRef = normErr*ones(size(error1));

rezu1 = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2 = [iterVec resRedVec, resPhyVec, error1, errorRef];
name = sprintf('output/historyCGNR_CHDG_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver CGNE']);
[resRedVec, resPhyVec, error2] = solverCGNEredu_DG(mesh, dofm, sysA, tol, iMax, iOut, @computeNormError2D_DG);

iterVec = (0:iOut:iMax)';
errorRef = normErr*ones(size(error2));

rezu1 = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2 = [iterVec resRedVec, resPhyVec, error2, errorRef];
name = sprintf('output/historyCGNE_CHDG_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver GMRES']);
[resRedVec, resPhyVec, error3] = solverGMRESredu_DG(mesh, dofm, sysA, tol, iMax, iOut, @computeNormError2D_DG);

iterVec = (0:iOut:iMax)';
errorRef = normErr*ones(size(error3));

rezu1 = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2 = [iterVec resRedVec, resPhyVec, error3, errorRef];
name = sprintf('output/historyGMRES_CHDG_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorNum = normErr*ones(size(error3));
errorProj = normProjErr*ones(size(error3));

figure;
hold off
semilogy(iterVec,error0,'DisplayName','Richardson');
hold on
semilogy(iterVec,error1,'DisplayName','CGNR');
semilogy(iterVec,error2,'DisplayName','CGNE');
semilogy(iterVec,error3,'DisplayName','GMRES');
semilogy(iterVec,errorNum,'k--','DisplayName','Numerical error');
semilogy(iterVec,errorProj,'k:','DisplayName','Projection error');
box on;
grid on;
legend('Location','southwest');
xlabel('Iteration');
ylabel('Relative error');
axis([0 iMax 0.0025 1]);

% figure;
% hold off
% semilogy(iterVec,resPhyVec,'-o','DisplayName','Relative residual (Phy)');
% hold on
% semilogy(iterVec,resRedVec,'-x','DisplayName','Relative residual (Red)');
% semilogy(iterVec,error    ,'-o','DisplayName','Relative L2-error');
% semilogy(iterVec,errorRef ,'k--','DisplayName','Relative L2-error (Ref)');
% box on;
% grid on;
% legend('Location','southwest');
% xlabel('Iteration');
% ylabel('Value');
% title(['CHDG ' benchmark ' — k=' num2str(k/pi) 'pi — h=' num2str(degree) ' — h=' num2str(h)]);
% axis([0 iMax 1e-10 1]);

end