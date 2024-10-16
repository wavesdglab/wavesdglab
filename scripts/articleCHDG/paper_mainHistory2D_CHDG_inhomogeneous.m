% close all;
clear all;

global h1 h2
global omega eta1 eta2 k1 k2 c1 c2 rho1 rho2
global rho c eta k

degree = 3;
tol = 1e-100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLANE WAVE
% iMax = 1000; iOut = 50; degree = 3;
% benchmark = 'open_heterogeneous';
% h1 = 1/16; h2 = 1/34; omega = 15*pi;
% rho1 = 1; c1 = 1; rho2 = 1; c2 = 1/2;
% eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
% run(benchmark,degree,tol,iMax,iOut);

% CAVITY 
iMax = 1000; iOut = 50; degree = 3;
benchmark = 'disk_heterogeneous';
h1 = 0.1; h2 = 0.075; omega = 10*pi;         % Case (1)
% h = 0.055; omega = 36.14;      % Case (2)
rho1 = 1; c1 = 1; rho2 = 1; c2 = 2/3;
eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
run(benchmark,degree,tol,iMax,iOut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,tol,iMax,iOut)
global h1 h2 k1 k2

% Build mesh and dofManager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method CHDG (' benchmark ')']);
disp(['---------------------------------------------------------']);
disp(['    h1                  ' num2str(h1)]);
disp(['    h2                  ' num2str(h2)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

BASIS = 0; PREC = 1;

[solA, sysA] = computeSolNum2D_CHDG_sym(mesh, dofm, PREC, 2);
[solB, sysB] = computeSolNum2D_CHDG_sym(mesh, dofm, PREC, 1);
[solC, sysC] = computeSolNum2D_CHDG_upw(mesh, dofm, PREC);
[solD, sysD] = computeSolNum2D_HDG_sym(mesh, dofm, BASIS, PREC);

normErrA = computeNormError2D_DG_ALL(mesh, dofm, solA)
normErrB = computeNormError2D_DG_ALL(mesh, dofm, solB)
normErrC = computeNormError2D_DG_ALL(mesh, dofm, solC)
normErrD = computeNormError2D_DG_ALL(mesh, dofm, solD)

solP = computeSolProjL2_2D_DG(mesh, dofm);
normProjErr = computeNormError2D_DG_ALL(mesh, dofm, solP);

% disp(['    L2-Error (numSol)   ' num2str(normErrA,'%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(normProjErr,'%1.2e')]);
disp(['---------------------------------------------------------']);

% writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
writeField2D(dofm, mesh, solC, 'output/solNum.pos', "solNum");
system('gmsh output/solNum.pos&');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver Richardson']);
alpha = 1;
[resRedVecA, resPhyVecA, error0A] = solverRichardsonRedu_DG(mesh, dofm, sysA, tol, iMax, iOut, alpha, @computeNormError2D_DG_ALL);
[resRedVecB, resPhyVecB, error0B] = solverRichardsonRedu_DG(mesh, dofm, sysB, tol, iMax, iOut, alpha, @computeNormError2D_DG_ALL);
[resRedVecC, resPhyVecC, error0C] = solverRichardsonRedu_DG(mesh, dofm, sysC, tol, iMax, iOut, alpha, @computeNormError2D_DG_ALL);

iterVec = (0:iOut:iMax)';
errorRefA = normErrA*ones(size(error0A));
errorRefB = normErrB*ones(size(error0B));
errorRefC = normErrC*ones(size(error0C));

rezu1A = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2A = [iterVec resRedVecA, resPhyVecA, error0A, errorRefA];
name = sprintf('output/historyRichardson_CHDG_2ndorder_%s_p%i_h1%g_h2%g_k1%g_k2%g.csv', benchmark, degree, h1, h2, k1, k2);
writematrix([rezu1A ; rezu2A], name, 'Delimiter', 'semi');

rezu1B = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2B = [iterVec resRedVecB, resPhyVecB, error0B, errorRefB];
name = sprintf('output/historyRichardson_CHDG_0thorder_%s_p%i_h1%g_h2%g_k1%g_k2%g.csv', benchmark, degree, h1, h2, k1, k2);
writematrix([rezu1B ; rezu2B], name, 'Delimiter', 'semi');

rezu1C = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2C = [iterVec resRedVecC, resPhyVecC, error0C, errorRefC];
name = sprintf('output/historyRichardson_CHDG_upwind_%s_p%i_h1%g_h2%g_k1%g_k2%g.csv', benchmark, degree, h1, h2, k1, k2);
writematrix([rezu1C ; rezu2C], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver CGNR']);
[resRedVecA, resPhyVecA, error1A] = solverCGNRredu_DG(mesh, dofm, sysA, tol, iMax, iOut, @computeNormError2D_DG_ALL);
[resRedVecB, resPhyVecB, error1B] = solverCGNRredu_DG(mesh, dofm, sysB, tol, iMax, iOut, @computeNormError2D_DG_ALL);
[resRedVecC, resPhyVecC, error1C] = solverCGNRredu_DG(mesh, dofm, sysC, tol, iMax, iOut, @computeNormError2D_DG_ALL);
[resRedVecD, resPhyVecD, error1D] = solverCGNRredu_DG(mesh, dofm, sysD, tol, iMax, iOut, @computeNormError2D_DG_ALL);

iterVec = (0:iOut:iMax)';
errorRefA = normErrA*ones(size(error1A));
errorRefB = normErrB*ones(size(error1B));
errorRefC = normErrC*ones(size(error1C));
errorRefD = normErrD*ones(size(error1D));

rezu1A = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2A = [iterVec resRedVecA, resPhyVecA, error1A, errorRefA];
name = sprintf('output/historyCGNR_CHDG_2ndorder_%s_p%i_h1%g_h2%g_k1%g_k2%g.csv', benchmark, degree, h1, h2, k1, k2);
writematrix([rezu1A ; rezu2A], name, 'Delimiter', 'semi');

rezu1B = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2B = [iterVec resRedVecB, resPhyVecB, error1B, errorRefB];
name = sprintf('output/historyCGNR_CHDG_0thorder_%s_p%i_h1%g_h2%g_k1%g_k2%g.csv', benchmark, degree, h1, h2, k1, k2);
writematrix([rezu1B ; rezu2B], name, 'Delimiter', 'semi');

rezu1C = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2C = [iterVec resRedVecC, resPhyVecC, error1C, errorRefC];
name = sprintf('output/historyCGNR_CHDG_upwind_%s_p%i_h1%g_h2%g_k1%g_k2%g.csv', benchmark, degree, h1, h2, k1, k2);
writematrix([rezu1C ; rezu2C], name, 'Delimiter', 'semi');

rezu1D = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2D = [iterVec resRedVecD, resPhyVecD, error1D, errorRefD];
name = sprintf('output/historyCGNR_HDG_highorder_%s_p%i_h1%g_h2%g_k1%g_k2%g.csv', benchmark, degree, h1, h2, k1, k2);
writematrix([rezu1D ; rezu2D], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver GMRES']);
[resRedVecA, resPhyVecA, error3A] = solverGMRESredu_DG(mesh, dofm, sysA, tol, iMax, iOut, @computeNormError2D_DG_ALL);
[resRedVecB, resPhyVecB, error3B] = solverGMRESredu_DG(mesh, dofm, sysB, tol, iMax, iOut, @computeNormError2D_DG_ALL);
[resRedVecC, resPhyVecC, error3C] = solverGMRESredu_DG(mesh, dofm, sysC, tol, iMax, iOut, @computeNormError2D_DG_ALL);
[resRedVecD, resPhyVecD, error3D] = solverGMRESredu_DG(mesh, dofm, sysD, tol, iMax, iOut, @computeNormError2D_DG_ALL);

iterVec = (0:iOut:iMax)';
errorRefA = normErrA*ones(size(error3A));
errorRefB = normErrB*ones(size(error3B));
errorRefC = normErrC*ones(size(error3C));
errorRefD = normErrD*ones(size(error3D));

rezu1A = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2A = [iterVec resRedVecA, resPhyVecA, error3A, errorRefA];
name = sprintf('output/historyGMRES_CHDG_2ndorder_%s_p%i_h1%g_h2%g_k1%g_k2%g.csv', benchmark, degree, h1, h2, k1, k2);
writematrix([rezu1A ; rezu2A], name, 'Delimiter', 'semi');

rezu1B = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2B = [iterVec resRedVecB, resPhyVecB, error3B, errorRefB];
name = sprintf('output/historyGMRES_CHDG_0thorder_%s_p%i_h1%g_h2%g_k1%g_k2%g.csv', benchmark, degree, h1, h2, k1, k2);
writematrix([rezu1B ; rezu2B], name, 'Delimiter', 'semi');

rezu1C = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2C = [iterVec resRedVecC, resPhyVecC, error3C, errorRefC];
name = sprintf('output/historyGMRES_CHDG_upwind_%s_p%i_h1%g_h2%g_k1%g_k2%g.csv', benchmark, degree, h1, h2, k1, k2);
writematrix([rezu1C ; rezu2C], name, 'Delimiter', 'semi');

rezu1D = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2D = [iterVec resRedVecD, resPhyVecD, error3D, errorRefD];
name = sprintf('output/historyGMRES_HDG_highorder_%s_p%i_h1%g_h2%g_k1%g_k2%g.csv', benchmark, degree, h1, h2, k1, k2);
writematrix([rezu1D ; rezu2D], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% errorNumA = normErrA*ones(size(error3A));
% errorProj = normProjErr*ones(size(error3A));

figure;
hold off
semilogy(iterVec,error0A,'-+y','DisplayName','CHDG (Sym-2) - Richardson');
hold on
semilogy(iterVec,error0B,'-xg','DisplayName','CHDG (Sym-0) - Richardson');
semilogy(iterVec,error0C,'-xb','DisplayName','CHDG (upwind) - Richardson');
semilogy(iterVec,error1B,'-oc','MarkerFaceColor','c','DisplayName','CHDG (Sym-0) - CGNR');
semilogy(iterVec,error3B,'-oc','MarkerFaceColor','w','DisplayName','CHDG (Sym-0) - GMRES');

semilogy(iterVec,error1A,'-oy','MarkerFaceColor','y','DisplayName','CHDG (Sym-2) - CGNR');
semilogy(iterVec,error3A,'-oy','MarkerFaceColor','w','DisplayName','CHDG (Sym-2) - GMRES');
semilogy(iterVec,error1C,'-ob','MarkerFaceColor','b','DisplayName','CHDG (upwind) - CGNR');
semilogy(iterVec,error3C,'-ob','MarkerFaceColor','w','DisplayName','CHDG (upwind) - GMRES');
semilogy(iterVec,error1D,'-or','MarkerFaceColor','r','DisplayName','HDG - CGNR');
hold on
semilogy(iterVec,error3D,'-or','MarkerFaceColor','w','DisplayName','HDG - GMRES');
% semilogy(iterVec,errorNumA,'k--','DisplayName','Numerical error');
% semilogy(iterVec,errorProj,'k:','DisplayName','Projection error');
box on;
grid on;
legend('Location','southwest');
xlabel('Iteration');
ylabel('Relative error');
axis([0 iMax 0.005 1]);

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