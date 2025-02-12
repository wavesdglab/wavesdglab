%close all;
clear;

global omega c rho h M theta phi v0

degree = 3;
tol = 1e-100;
iMax = 1000;
iOut = 50;
PREC = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Benchmark 'Convected + Cavity'
benchmark = 'open_convected';
omega = 15*pi; c = 1; rho = 1; h = 1/30;
M = 0.5; theta = 5*pi/4; phi = pi/4; v0 = [M*c*cos(theta), M*c*sin(theta)];
run(M,benchmark,degree,PREC,tol,iMax,iOut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(M,benchmark,degree,PREC,tol,iMax,iOut)

global c rho phi theta

% Build mesh and dofManager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Benchmark ' benchmark '']);
disp(['---------------------------------------------------------']);
disp(['    c                   ' num2str(c)]);
disp(['    rho                 ' num2str(rho)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

% [solDG, sysDG] = computeSolNum2D_DG_convected(mesh, dofm, PREC);
% [solCG, sysCG] = computeSolNum2D_CHDG_convected_v5(mesh, dofm, PREC);
% [solDG, sysDG] = computeSolNum2D_CHDG_convected_v4(mesh, dofm, PREC);
[solHDG, sysHDG] = computeSolNum2D_CHDG_convected_v6(mesh, dofm, PREC);
[solCHDG, sysCHDG] = computeSolNum2D_CHDG_convected_v2(mesh, dofm, PREC);
% normErrCG = computeNormError2D_DG_convected(mesh, dofm, solCG);
% normErrDG = computeNormError2D_DG_convected(mesh, dofm, solDG);
normErrCHDG = computeNormError2D_DG_convected(mesh, dofm, solCHDG);
normErrHDG = computeNormError2D_DG_convected(mesh, dofm, solHDG);

% disp(['    L2-Error (CG)       ' num2str(normErrCG,'%1.2e')]);
% disp(['    L2-Error (DG)       ' num2str(normErrDG,'%1.2e')]);
disp(['    L2-Error (CHDG)     ' num2str(normErrCHDG,'%1.2e')]);
disp(['    L2-Error (HDG)      ' num2str(normErrHDG,'%1.2e')]);

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

disp('    (CHDG)');
[resRedVecCHDG,  resPhyVecCHDG,  errorRichCHDG ] = solverRichardsonRedu_DG(mesh, dofm, sysCHDG,  tol, iMax, iOut, alpha, @computeNormError2D_DG_convected, []);
rezu = [iterVec, resRedVecCHDG, resPhyVecCHDG, errorRichCHDG, normErrCHDG*ones(size(iterVec))];
name = sprintf('output/historyRichardson_CHDG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (HDG)');
[resRedVecHDG,  resPhyVecHDG,  errorRichHDG ] = solverRichardsonRedu_DG(mesh, dofm, sysHDG,  tol, iMax, iOut, alpha, @computeNormError2D_DG_convected, []);
rezu = [iterVec, resRedVecHDG, resPhyVecHDG, errorRichHDG, normErrHDG*ones(size(iterVec))];
name = sprintf('output/historyRichardson_HDG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');
 
% % % disp('    (DG)');
% % % [resVecDG, resPhyVecDG, errorRichDG] = solverRichardsonRedu_DG(mesh, dofm, sysDG, tol, iMax, iOut, alpha, @computeNormError2D_DG_convected, []);
% % % rezu = [iterVec, resVecDG, resPhyVecDG, errorRichDG, normErrDG*ones(size(iterVec))];
% % % name = sprintf('output/historyRichardson_DG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
% % % writematrix([labels ; rezu], name, 'Delimiter', 'semi');

% disp('    (DG)');
% [resVecDG, resPhyVecDG, errorRichDG] = solverRichardson_DG(mesh, dofm, sysDG, tol, iMax, iOut, alpha, @computeNormError2D_DG_convected);
% rezu = [iterVec, resVecDG, resPhyVecDG, errorRichDG, normErrDG*ones(size(iterVec))];
% name = sprintf('output/historyRichardson_DG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
% writematrix([labels ; rezu], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--- Solver CGNR');

disp('    (CHDG)');
[resRedVecCHDG, resPhyVecCHDG, errorCgnrCHDG] = solverCGNRredu_DG(mesh, dofm, sysCHDG, tol, iMax, iOut, @computeNormError2D_DG_convected, []);
rezu = [iterVec, resRedVecCHDG, resPhyVecCHDG, errorCgnrCHDG, normErrCHDG*ones(size(iterVec))];
name = sprintf('output/historyCGNR_CHDG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (HDG)');
[resRedVecHDG, resPhyVecHDG, errorCgnrHDG] = solverCGNRredu_DG(mesh, dofm, sysHDG, tol, iMax, iOut, @computeNormError2D_DG_convected, []);
rezu = [iterVec, resRedVecHDG, resPhyVecHDG, errorCgnrHDG, normErrHDG*ones(size(iterVec))];
name = sprintf('output/historyCGNR_HDG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

% disp('    (DG)');
% [resRedVecDG, resPhyVecDG, errorCgnrDG] = solverCGNRredu_DG(mesh, dofm, sysDG, tol, iMax, iOut, @computeNormError2D_DG_convected, []);
% rezu = [iterVec, resRedVecDG, resPhyVecDG, errorCgnrDG, normErrDG*ones(size(iterVec))];
% name = sprintf('output/historyCGNR_DG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
% writematrix([labels ; rezu], name, 'Delimiter', 'semi');
% 
% disp('    (CG)');
% [resRedVecCG, resPhyVecCG, errorCgnrCG] = solverCGNRredu_DG(mesh, dofm, sysCG, tol, iMax, iOut, @computeNormError2D_DG_convected, []);
% rezu = [iterVec, resRedVecCG, resPhyVecCG, errorCgnrCG, normErrDG*ones(size(iterVec))];
% name = sprintf('output/historyCGNR_CG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
% writematrix([labels ; rezu], name, 'Delimiter', 'semi');

% disp('    (DG)');
% [resVecDG, errorCgnrDG] = solverCGNR(mesh, dofm, sysDG, tol, iMax, iOut, @computeNormError2D_DG_convected);
% resPhyVecDG = resVecDG;
% rezu = [iterVec, resVecDG, resPhyVecDG, errorCgnrDG, normErrDG*ones(size(iterVec))];
% name = sprintf('output/historyCGNR_DG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
% writematrix([labels ; rezu], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--- Solver GMRES');

disp('    (CHDG)');
[resRedVecCHDG, resPhyVecCHDG, errorGmresCHDG] = solverGMRESredu_DG(mesh, dofm, sysCHDG, tol, iMax, iOut, @computeNormError2D_DG_convected, []);
rezu = [iterVec, resRedVecCHDG, resPhyVecCHDG, errorGmresCHDG, normErrCHDG*ones(size(iterVec))];
name = sprintf('output/historyGMRES_CHDG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

disp('    (HDG)');
[resRedVecHDG, resPhyVecHDG, errorGmresHDG] = solverGMRESredu_DG(mesh, dofm, sysHDG, tol, iMax, iOut, @computeNormError2D_DG_convected, []);
rezu = [iterVec, resRedVecHDG, resPhyVecHDG, errorGmresHDG, normErrCHDG*ones(size(iterVec))];
name = sprintf('output/historyGMRES_HDG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
writematrix([labels ; rezu], name, 'Delimiter', 'semi');

% % % disp('    (DG)');
% % % [resRedVecDG, resPhyVecDG, errorGmresDG] = solverGMRESredu_DG(mesh, dofm, sysDG, tol, iMax, iOut, @computeNormError2D_DG_convected, []);
% % % rezu = [iterVec, resRedVecDG, resPhyVecDG, errorGmresDG, normErrDG*ones(size(iterVec))];
% % % name = sprintf('output/historyGMRES_DG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
% % % writematrix([labels ; rezu], name, 'Delimiter', 'semi');

% disp('    (DG)');
% [resVecDG, errorGmresDG] = solverGMRES(mesh, dofm, sysDG, tol, iMax, iOut, @computeNormError2D_DG_convected);
% resPhyVecDG = resVecDG;
% rezu = [iterVec, resVecDG, resPhyVecDG, errorGmresDG, normErrDG*ones(size(iterVec))];
% name = sprintf('output/historyGMRES_DG_%s_p%i_M%g_prec%i_c_%g_rho_%g_phi_%g_theta_%g.csv', benchmark, degree, M, PREC, c, rho, phi, theta);
% writematrix([labels ; rezu], name, 'Delimiter', 'semi');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();
hold off
% semilogy(iterVec,errorRichDG, '-xg','DisplayName','DG - Richardson'); hold on; pause(0.1);
semilogy(iterVec,errorRichCHDG, '-xr','DisplayName','CHDG - Richardson'); hold on; pause(0.1);
semilogy(iterVec,errorRichHDG, '-xb','DisplayName','HDG - Richardson'); hold on; pause(0.1);
% semilogy(iterVec,errorCgnrCG,   '-oc','MarkerFaceColor','w','DisplayName','CG - CGNR'); hold on; pause(0.1);
% semilogy(iterVec,errorCgnrDG,   '-og','MarkerFaceColor','w','DisplayName','DG - CGNR'); hold on; pause(0.1);
semilogy(iterVec,errorCgnrCHDG,   '-or','MarkerFaceColor','w','DisplayName','CHDG - CGNR'); hold on; pause(0.1);
semilogy(iterVec,errorCgnrHDG,   '-ob','MarkerFaceColor','w','DisplayName','HDG - CGNR'); hold on; pause(0.1);
% semilogy(iterVec,errorGmresDG,   '-og','MarkerFaceColor','g','DisplayName','DG - GMRES'); hold on; pause(0.1);
semilogy(iterVec,errorGmresCHDG,   '-or','MarkerFaceColor','r','DisplayName','CHDG - GMRES'); hold on; pause(0.1);
semilogy(iterVec,errorGmresHDG,   '-ob','MarkerFaceColor','b','DisplayName','HDG - GMRES'); hold on; pause(0.1);
errorCHDG = normErrCHDG*ones(size(iterVec));
% errorProj = normErrProj*ones(size(iterVec));
semilogy(iterVec,errorCHDG,'k--','DisplayName','Numerical error');
% semilogy(iterVec,errorProj,'k:','DisplayName','Projection error');

box on;
grid on;
legend('Location','southwest');
xlabel('Iteration');
ylabel('Relative error');
axis([0 iMax 0.0005 1]);

end