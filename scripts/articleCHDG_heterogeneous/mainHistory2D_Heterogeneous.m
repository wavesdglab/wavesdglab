%close all;
clear;

global omega c1 c2 rho1 rho2 h1 h2

degree = 3;
tol = 1e-100;
iMax = 1000;
iOut = 50;
PREC = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Benchmark 'Homogeneous + Open'
benchmark = 'open_heterogeneous';
omega = 15*pi; c1 = 1; c2 = 1; rho1 = 1; rho2 = 1; h1 = 1/16; h2 = h1;
run(benchmark,degree,PREC,tol,iMax,iOut);
omega = 30*pi; c1 = 1; c2 = 1; rho1 = 1; rho2 = 1; h1 = 1/34; h2 = h1;
run(benchmark,degree,PREC,tol,iMax,iOut);

% Benchmark 'Homogeneous + Cavity'
benchmark = 'disk_heterogeneous';
omega = 16.5;  c1 = 1; c2 = 1; rho1 = 1; rho2 = 1; h1 = 0.04; h2 = h1;
run(benchmark,degree,PREC,tol,iMax,iOut);
omega = 17;    c1 = 1; c2 = 1; rho1 = 1; rho2 = 1; h1 = 0.025; h2 = h1;
run(benchmark,degree,PREC,tol,iMax,iOut);

% Benchmark 'Heterogeneous + Open'
benchmark = 'open_heterogeneous';
omega = 15*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 2; h1 = 1/16; h2 = 1/34;
run(benchmark,degree,PREC,tol,iMax,iOut);
omega = 15*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 1; h1 = 1/16; h2 = 1/34;
run(benchmark,degree,PREC,tol,iMax,iOut);

% Benchmark 'Heterogeneous + Cavity'
benchmark = 'disk_heterogeneous';
omega = 10*pi; c1 = 1; c2 = 2/3; rho1 = 1; rho2 = 3/2; h1 = 1/12; h2 = 1/16;
run(benchmark,degree,PREC,tol,iMax,iOut);
omega = 10*pi; c1 = 1; c2 = 2/3; rho1 = 1; rho2 = 1; h1 = 1/12; h2 = 1/16;
run(benchmark,degree,PREC,tol,iMax,iOut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,PREC,tol,iMax,iOut
)
global k1 k2 eta1 eta2

% Build mesh and dofManager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Benchmark ' benchmark '']);
disp(['---------------------------------------------------------']);
disp(['    k1                  ' num2str(k1)]);
disp(['    k2                  ' num2str(k2)]);
disp(['    eta1                ' num2str(eta1)]);
disp(['    eta2                ' num2str(eta2)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

[solUpw, sysUpw] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, 'UPW', PREC);
normErrUpw = computeNormError2D_DG_heterogeneous(mesh, dofm, solUpw);
disp(['    L2-Error (CHDG UPW)   ' num2str(normErrUpw,'%1.2e')]);

[solSym, sysSym] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, 'SYM', PREC);
normErrSym = computeNormError2D_DG_heterogeneous(mesh, dofm, solSym);
disp(['    L2-Error (CHDG SYM)   ' num2str(normErrSym,'%1.2e')]);

[solSym2, sysSym2] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, 'SYM2', PREC);
normErrSym2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solSym2);
disp(['    L2-Error (CHDG SYM2)  ' num2str(normErrSym2,'%1.2e')]);

[solHdgUpw, sysHdgUpw] = computeSolNum2D_HDG_heterogeneous(mesh, dofm, 'UPW');
normErrHdgUpw = computeNormError2D_DG_heterogeneous(mesh, dofm, solHdgUpw);
disp(['    L2-Error (HDG UPW)    ' num2str(normErrHdgUpw,'%1.2e')]);

[solHdgSym, sysHdgSym] = computeSolNum2D_HDG_heterogeneous(mesh, dofm, 'SYM');
normErrHdgSym = computeNormError2D_DG_heterogeneous(mesh, dofm, solHdgSym);
disp(['    L2-Error (HDG SYM)    ' num2str(normErrHdgSym,'%1.2e')]);

solP = computeSolProjL2_2D_DG(mesh, dofm);
normErrProj = computeNormError2D_DG_heterogeneous(mesh, dofm, solP);
disp(['    L2-Error (projSol)  ' num2str(normErrProj,'%1.2e')]);

disp('---------------------------------------------------------');

% writeField2D(dofm, mesh, solSym2, 'output/solNum.pos', "solNum");
% writeField2D(dofm, mesh, solUpw, 'output/solNum.pos', "solNum");
% system('gmsh output/solNum.pos&');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iterVec = (0:iOut:iMax)';
labels = ["iter" "resRed" "resPhy" "error" "errorRef"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--- Solver Richardson');

alpha = 1;

disp('    (CHDG UPW)');
[resRedVecUpw,  resPhyVecUpw,  errorRichUpw ] = solverRichardsonRedu_DG(mesh, dofm, sysUpw,  tol, iMax, iOut, alpha, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecUpw, resPhyVecUpw, errorRichUpw, normErrUpw*ones(size(iterVec))];
name = sprintf('output/historyRichardson_CHDG_UPW_%s_p%i_prec%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, PREC, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (CHDG SYM)');
[resRedVecSym, resPhyVecSym, errorRichSym] = solverRichardsonRedu_DG(mesh, dofm, sysSym,  tol, iMax, iOut, alpha, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecSym, resPhyVecSym, errorRichSym, normErrSym*ones(size(iterVec))];
name = sprintf('output/historyRichardson_CHDG_SYM_%s_p%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (CHDG SYM2)');
[resRedVecSym2, resPhyVecSym2, errorRichSym2] = solverRichardsonRedu_DG(mesh, dofm, sysSym2, tol, iMax, iOut, alpha, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecSym2, resPhyVecSym2, errorRichSym2, normErrSym2*ones(size(iterVec))];
name = sprintf('output/historyRichardson_CHDG_SYM2_%s_p%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--- Solver CGNR');

disp('    (CHDG UPW)');
[resRedVecUpw, resPhyVecUpw, errorCgnrUpw] = solverCGNRredu_DG(mesh, dofm, sysUpw, tol, iMax, iOut, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecUpw, resPhyVecUpw, errorCgnrUpw, normErrUpw*ones(size(iterVec))];
name = sprintf('output/historyCGNR_CHDG_UPW_%s_p%i_prec%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, PREC, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (CHDG SYM)');
[resRedVecSym, resPhyVecSym, errorCgnrSym] = solverCGNRredu_DG(mesh, dofm, sysSym, tol, iMax, iOut, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecSym, resPhyVecSym, errorCgnrSym, normErrSym*ones(size(iterVec))];
name = sprintf('output/historyCGNR_CHDG_SYM_%s_p%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (CHDG SYM2)');
[resRedVecSym2, resPhyVecSym2, errorCgnrSym2] = solverCGNRredu_DG(mesh, dofm, sysSym2, tol, iMax, iOut, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecSym2, resPhyVecSym2, errorCgnrSym2, normErrSym2*ones(size(iterVec))];
name = sprintf('output/historyCGNR_CHDG_SYM2_%s_p%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (HDG UPW)');
[resRedVecHdgUpw, resPhyVecHdgUpw, errorCgnrHdgUpw] = solverCGNRredu_DG(mesh, dofm, sysHdgUpw, tol, iMax, iOut, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecHdgUpw, resPhyVecHdgUpw, errorCgnrHdgUpw, normErrHdgSym*ones(size(iterVec))];
name = sprintf('output/historyCGNR_HDG_UPW_%s_p%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (HDG SYM)');
[resRedVecHdgSym, resPhyVecHdgSym, errorCgnrHdgSym] = solverCGNRredu_DG(mesh, dofm, sysHdgSym, tol, iMax, iOut, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecHdgSym, resPhyVecHdgSym, errorCgnrHdgSym, normErrHdgSym*ones(size(iterVec))];
name = sprintf('output/historyCGNR_HDG_SYM_%s_p%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--- Solver GMRES');

disp('    (CHDG UPW)');
[resRedVecUpw, resPhyVecUpw, errorGmresUpw] = solverGMRESredu_DG(mesh, dofm, sysUpw, tol, iMax, iOut, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecUpw, resPhyVecUpw, errorGmresUpw, normErrUpw*ones(size(iterVec))];
name = sprintf('output/historyGMRES_CHDG_UPW_%s_p%i_prec%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, PREC, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (CHDG SYM)');
[resRedVecSym, resPhyVecSym, errorGmresSym] = solverGMRESredu_DG(mesh, dofm, sysSym, tol, iMax, iOut, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecSym, resPhyVecSym, errorGmresSym, normErrSym*ones(size(iterVec))];
name = sprintf('output/historyGMRES_CHDG_SYM_%s_p%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (CHDG SYM2)');
[resRedVecSym2, resPhyVecSym2, errorGmresSym2] = solverGMRESredu_DG(mesh, dofm, sysSym2, tol, iMax, iOut, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecSym2, resPhyVecSym2, errorGmresSym2, normErrSym2*ones(size(iterVec))];
name = sprintf('output/historyGMRES_CHDG_SYM2_%s_p%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (HDG UPW)');
[resRedVecHdgUpw, resPhyVecHdgUpw, errorGmresHdgUpw] = solverGMRESredu_DG(mesh, dofm, sysHdgUpw, tol, iMax, iOut, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecHdgUpw, resPhyVecHdgUpw, errorGmresHdgUpw, normErrHdgSym*ones(size(iterVec))];
name = sprintf('output/historyGMRES_HDG_UPW_%s_p%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (HDG SYM)');
[resRedVecHdgSym, resPhyVecHdgSym, errorGmresHdgSym] = solverGMRESredu_DG(mesh, dofm, sysHdgSym, tol, iMax, iOut, @computeNormError2D_DG_heterogeneous, []);
rezu = [iterVec, resRedVecHdgSym, resPhyVecHdgSym, errorGmresHdgSym, normErrHdgSym*ones(size(iterVec))];
name = sprintf('output/historyGMRES_HDG_SYM_%s_p%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, k1, k2, eta1, eta2);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(1);
% hold off
% semilogy(iterVec,errorRichUpw, '-xb','DisplayName','CHDG (UPW) - Richardson'); hold on; pause(0.1);
% semilogy(iterVec,errorRichSym, '-xk','DisplayName','CHDG (SYM) - Richardson'); hold on; pause(0.1);
% semilogy(iterVec,errorRichSym2,'-xm','DisplayName','CHDG (SYM2) - Richardson'); hold on; pause(0.1);
% semilogy(iterVec,errorCgnrUpw,   '-ob','MarkerFaceColor','b','DisplayName','CHDG (UPW) - CGNR'); hold on; pause(0.1);
% semilogy(iterVec,errorCgnrSym,   '-ok','MarkerFaceColor','k','DisplayName','CHDG (SYM) - CGNR'); hold on; pause(0.1);
% semilogy(iterVec,errorCgnrSym2,  '-om','MarkerFaceColor','m','DisplayName','CHDG (SYM2) - CGNR'); hold on; pause(0.1);
% semilogy(iterVec,errorCgnrHdgUpw,'-oc','MarkerFaceColor','c','DisplayName','HDG (UPW) - CGNR'); hold on; pause(0.1);
% semilogy(iterVec,errorCgnrHdgSym,'-or','MarkerFaceColor','r','DisplayName','HDG (SYM) - CGNR'); hold on; pause(0.1);
% semilogy(iterVec,errorGmresUpw,   '-ob','MarkerFaceColor','w','DisplayName','CHDG (UPW) - GMRES'); hold on; pause(0.1);
% semilogy(iterVec,errorGmresSym,   '-ok','MarkerFaceColor','w','DisplayName','CHDG (SYM) - GMRES'); hold on; pause(0.1);
% semilogy(iterVec,errorGmresSym2,  '-om','MarkerFaceColor','w','DisplayName','CHDG (SYM2) - GMRES'); hold on; pause(0.1);
% semilogy(iterVec,errorGmresHdgUpw,'-oc','MarkerFaceColor','w','DisplayName','HDG (UPW) - GMRES'); hold on; pause(0.1);
% semilogy(iterVec,errorGmresHdgSym,'-or','MarkerFaceColor','w','DisplayName','HDG (SYM) - GMRES'); hold on; pause(0.1);
% 
% errorSym = normErrSym*ones(size(iterVec));
% errorProj = normErrProj*ones(size(iterVec));
% semilogy(iterVec,errorSym,'k--','DisplayName','Numerical error (SYM)');
% semilogy(iterVec,errorProj,'k:','DisplayName','Projection error');
% 
% box on;
% grid on;
% legend('Location','southwest');
% xlabel('Iteration');
% ylabel('Relative error');
% axis([0 iMax 0.005 1]);

end