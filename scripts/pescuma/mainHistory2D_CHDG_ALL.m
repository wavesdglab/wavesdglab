% close all;
clear all;

global h
global omega eta1 eta2 k1 k2 c1 c2 rho1 rho2
global rho c eta k

degree = 3;
A = 2;              % order of numerical fluxes
B = 2;              % order of transmission variables
tol = 1e-100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLANE WAVE IN HOMOGENEOUS MEDIUM WITH EXACT DIRICHLET B.C.
% iMax = 1000; iOut = 50;
% benchmark = 'open_heterogeneous';
% 
% degree = 3; h = 1/10; omega = 15*pi;
% rho1 = 1; c1 = 1.3; rho2 = 1; c2 = 0.8;
% eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
% run(benchmark,degree,BASIS,PREC,A,B,tol,iMax,iOut);
% 
% degree = 3; h = 1/10; omega = 15*pi;
% rho1 = 1; c1 = 1.3; rho2 = 1; c2 = 0.8;
% eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
% run(benchmark,degree,BASIS,PREC,A,B,tol,iMax,iOut);

% PLANE WAVE IN INHOMOGENEOUS MEDIUM WITH EXACT ROBIN B.C.
iMax = 1000; iOut = 50;
benchmark = 'open_heterogeneous';

degree = 3; h = 1/34; omega = 15*pi;
rho1 = 1; c1 = 1; rho2 = 1; c2 = 0.5;
eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
run(benchmark,degree,A,B,tol,iMax,iOut);

% degree = 3; h = 1/10; omega = 15*pi;
% rho1 = 1; c1 = 1.3; rho2 = 1; c2 = 0.8;
% eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
% run(benchmark,degree,BASIS,PREC,A,B,tol,iMax,iOut);

% BENCH WAVEGUIDE
% iMax = 4000; iOut = 200;
% benchmark = 'waveguide'; k = 6*pi; h = 1/8;
% run(benchmark,degree,tau,BASIS,PREC,tol,iMax,iOut);
% benchmark = 'waveguide'; k = 12*pi; h = 1/17;
% run(benchmark,degree,tau,BASIS,PREC,tol,iMax,iOut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,A,B,tol,iMax,iOut)
global h
global omega eta1 eta2 k1 k2 c1 c2 rho1 rho2
global rho c eta k

% Build mesh and dofManager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);
setParameters(mesh,benchmark);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method CHDG (' benchmark ')']);
disp(['---------------------------------------------------------']);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

BASIS = 0; PREC = 1;

[solA, sysA] = computeSolNum2D_CHDG_ALL(mesh, dofm, PREC, A, B);

[~, sysB] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, BASIS, PREC);

[~, sysC] = computeSolNum2D_HDG_ALL(mesh, dofm, BASIS, PREC);

% [~, sysD] = computeSolNum2D_HDG_heterogeneous(mesh, dofm, BASIS, PREC);

normErr= computeNormError2D_DG_ALL(mesh, dofm, solA);

solP = computeSolProjL2_2D_DG(mesh, dofm);
normProjErr = computeNormError2D_DG_ALL(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(normErr,'%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(normProjErr,'%1.2e')]);
disp(['---------------------------------------------------------']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver Richardson']);
alpha = 1;
[resRedVecA, resPhyVecA, error0A] = solverRichardsonRedu_DG(mesh, dofm, sysA, tol, iMax, iOut, alpha, @computeNormError2D_DG_ALL);
[resRedVecB, resPhyVecB, error0B] = solverRichardsonRedu_DG(mesh, dofm, sysB, tol, iMax, iOut, alpha, @computeNormError2D_DG_ALL);

iterVec = (0:iOut:iMax)';
errorRefA = normErr*ones(size(error0A));
errorRefB = normErr*ones(size(error0B));

rezu1A = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2A = [iterVec resRedVecA, resPhyVecA, error0A, errorRefA];
name = sprintf('output/historyRichardson_CHDG_highorder_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
writematrix([rezu1A ; rezu2A], name, 'Delimiter', 'semi');

rezu1B = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2B = [iterVec resRedVecB, resPhyVecB, error0B, errorRefB];
name = sprintf('output/historyRichardson_CHDG_upwind_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
writematrix([rezu1B ; rezu2B], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver CGNR']);
[resRedVecA, resPhyVecA, error1A] = solverCGNRredu_DG(mesh, dofm, sysA, tol, iMax, iOut, @computeNormError2D_DG_ALL);
[resRedVecB, resPhyVecB, error1B] = solverCGNRredu_DG(mesh, dofm, sysB, tol, iMax, iOut, @computeNormError2D_DG_ALL);
[resRedVecC, resPhyVecC, error1C] = solverCGNRredu_DG(mesh, dofm, sysC, tol, iMax, iOut, @computeNormError2D_DG_ALL);

iterVec = (0:iOut:iMax)';
errorRefA = normErr*ones(size(error1A));
errorRefB = normErr*ones(size(error1B));
errorRefC = normErr*ones(size(error1C));

rezu1A = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2A = [iterVec resRedVecA, resPhyVecA, error1A, errorRefA];
name = sprintf('output/historyCGNR_CHDG_highorder_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
writematrix([rezu1A ; rezu2A], name, 'Delimiter', 'semi');

rezu1B = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2B = [iterVec resRedVecB, resPhyVecB, error1B, errorRefB];
name = sprintf('output/historyCGNR_CHDG_upwind_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
writematrix([rezu1B ; rezu2B], name, 'Delimiter', 'semi');

rezu1C = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2C = [iterVec resRedVecC, resPhyVecC, error1C, errorRefC];
name = sprintf('output/historyCGNR_HDG_highorder_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
writematrix([rezu1C ; rezu2C], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(['--- Solver CGNE']);
% [resRedVec, resPhyVec, error2] = solverCGNEredu_DG(mesh, dofm, sysA, tol, iMax, iOut, @computeNormError2D_DG);
% 
% iterVec = (0:iOut:iMax)';
% errorRef = normErr*ones(size(error2));
% 
% rezu1 = ["iter" "resRed" "resPhy" "error" "errorRef"];
% rezu2 = [iterVec resRedVec, resPhyVec, error2, errorRef];
% name = sprintf('output/historyCGNE_CHDG_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
% writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver GMRES']);
[resRedVecA, resPhyVecA, error3A] = solverGMRESredu_DG(mesh, dofm, sysA, tol, iMax, iOut, @computeNormError2D_DG_ALL);
[resRedVecB, resPhyVecB, error3B] = solverGMRESredu_DG(mesh, dofm, sysB, tol, iMax, iOut, @computeNormError2D_DG_ALL);
[resRedVecC, resPhyVecC, error3C] = solverGMRESredu_DG(mesh, dofm, sysC, tol, iMax, iOut, @computeNormError2D_DG_ALL);
% [resRedVecD, resPhyVecD, error3D] = solverGMRESredu_DG(mesh, dofm, sysD, tol, iMax, iOut, @computeNormError2D_DG_ALL);

iterVec = (0:iOut:iMax)';
errorRefA = normErr*ones(size(error3A));
errorRefB = normErr*ones(size(error3B));
errorRefC = normErr*ones(size(error3C));
% errorRefD = normErr*ones(size(error3D));

rezu1A = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2A = [iterVec resRedVecA, resPhyVecA, error3A, errorRefA];
name = sprintf('output/historyGMRES_CHDG_highorder_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
writematrix([rezu1A ; rezu2A], name, 'Delimiter', 'semi');

rezu1B = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2B = [iterVec resRedVecB, resPhyVecB, error3B, errorRefB];
name = sprintf('output/historyGMRES_CHDG_upwind_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
writematrix([rezu1B ; rezu2B], name, 'Delimiter', 'semi');

rezu1C = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2C = [iterVec resRedVecC, resPhyVecC, error3C, errorRefC];
name = sprintf('output/historyGMRES_HDG_highorder_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
writematrix([rezu1C ; rezu2C], name, 'Delimiter', 'semi');

% rezu1D = ["iter" "resRed" "resPhy" "error" "errorRef"];
% rezu2D = [iterVec resRedVecD, resPhyVecD, error3D, errorRefD];
% name = sprintf('output/historyGMRES_HDG_upwind_%s_p%i_h%g_k1%g_k2%g.csv', benchmark, degree, h, k1, k2);
% writematrix([rezu1D ; rezu2D], name, 'Delimiter', 'semi');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorNum = normErr*ones(size(error3A));
errorProj = normProjErr*ones(size(error3A));

figure;
hold off
semilogy(iterVec,error0A,'-*g','DisplayName','CHDG (high order) - Richardson');
hold on
semilogy(iterVec,error0B,'-xg','DisplayName','CHDG (upwind) - Richardson');
semilogy(iterVec,error1B,'-oc','MarkerFaceColor','c','DisplayName','CHDG (upwind) - CGNR');
semilogy(iterVec,error3B,'-oc','MarkerFaceColor','w','DisplayName','CHDG (upwind) - GMRES');
semilogy(iterVec,error1A,'-ob','MarkerFaceColor','b','DisplayName','CHDG (high order) - CGNR');
semilogy(iterVec,error3A,'-ob','MarkerFaceColor','w','DisplayName','CHDG (high order) - GMRES');
semilogy(iterVec,error1C,'-or','MarkerFaceColor','r','DisplayName','HDG - CGNR');
semilogy(iterVec,error3C,'-or','MarkerFaceColor','w','DisplayName','HDG - GMRES');
semilogy(iterVec,errorNum,'k--','DisplayName','Numerical error');
% semilogy(iterVec,errorProj,'k:','DisplayName','Projection error');
box on;
grid on;
legend('Location','southwest');
xlabel('Iteration');
ylabel('Relative error');
axis([0 iMax 0.05 1]);

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