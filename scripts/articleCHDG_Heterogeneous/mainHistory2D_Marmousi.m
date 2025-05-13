% close;
clear;

format long e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Benchmark
benchmark = 'geophysics_marmousi';

% freq = 5;
% iMax = 100;
% iOut = 10;
% iRestart = 10;
% tol = 1e-100;

freq = 30;
iMax = 10000;
iOut = 200;
iRestart = 10;
tol = 1e-100;

iterVec = (0:iOut:iMax)';
labels = ["iter" "resRed" "resPhy" "error"];

% Parameters
global omega nLambda
omega = 2*pi*freq;
degree = 3;
nLambda = 10/(degree+1);

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

METHOD = 'HDG';
FLUX = 'SYM';
PREC = 1;

if strcmp(METHOD,'CHDG')
    [sol, sys] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, FLUX, PREC);
end
if strcmp(METHOD,'HDG')
    [sol, sys] = computeSolNum2D_HDG_heterogeneous(mesh, dofm, FLUX);
end
xRef = sol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

valfig = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--- Solver Jacobi');
[resRedVec, resPhyVec, error] = solverRichardsonRedu_DG(mesh, dofm, sys, tol, iMax, iOut, 1, @computeNormError2D_DG_heterogeneous, sol);
rezuRich = [iterVec resRedVec, resPhyVec, error];
name = sprintf('output/historyJacobi_%s_%s_Marmousi_p%i_prec%i_omega%g.csv', METHOD, FLUX, degree, PREC, omega);
writematrix([labels ; rezuRich], name, 'Delimiter', 'semi');

figure(valfig+1); hold off; semilogy(iterVec,resRedVec,'DisplayName','Jacobi'); hold on; legend();
figure(valfig+2); hold off; semilogy(iterVec,resPhyVec,'DisplayName','Jacobi'); hold on; legend();
figure(valfig+3); hold off; semilogy(iterVec,error,'DisplayName','Jacobi'); hold on; legend();
pause(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--- Solver CGNR');
[resRedVec, resPhyVec, error] = solverCGNRredu_DG(mesh, dofm, sys, tol, iMax, iOut, @computeNormError2D_DG_heterogeneous, sol);
rezuCGNR = [iterVec resRedVec, resPhyVec, error];
name = sprintf('output/historyCGNR_%s_%s_Marmousi_p%i_prec%i_omega%g.csv', METHOD, FLUX, degree, PREC, omega);
writematrix([labels ; rezuCGNR], name, 'Delimiter', 'semi');

figure(valfig+1); hold off; semilogy(iterVec,resRedVec,'DisplayName','CGNR'); hold on; legend();
figure(valfig+2); hold off; semilogy(iterVec,resPhyVec,'DisplayName','CGNR'); hold on; legend();
figure(valfig+3); hold off; semilogy(iterVec,error,'DisplayName','CGNR'); hold on; legend();
pause(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--- Solver GMRES');
[resRedVec, resPhyVec, error] = solverGMRESredu_DG_restart(mesh, dofm, sys, tol, iRestart, iMax, iOut, @computeNormError2D_DG_heterogeneous, sol);
rezuGMRES = [iterVec resRedVec, resPhyVec, error];
name = sprintf('output/historyGMRES_%s_%s_Marmousi_p%i_prec%i_omega%g_restart%i.csv', METHOD, FLUX, degree, PREC, omega, iRestart);
writematrix([labels ; rezuGMRES], name, 'Delimiter', 'semi');

figure(valfig+1); hold off; semilogy(iterVec,resRedVec,'DisplayName','GMRES'); hold on; legend();
figure(valfig+2); hold off; semilogy(iterVec,resPhyVec,'DisplayName','GMRES'); hold on; legend();
figure(valfig+3); hold off; semilogy(iterVec,error,'DisplayName','GMRES'); hold on; legend();
pause(1);