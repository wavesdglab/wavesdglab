%close all;
clear;

global k h

PREC = 0;
tol = 1e-100;

% BENCH FREE SPACE
iMax = 300; iOut = 10;
benchmark = 'open'; k = 15*pi; h = 1/8; degree = 5;
run(benchmark,degree,PREC,tol,iMax,iOut);
% BENCH CAVITY
iMax = 600; iOut = 20;
benchmark = 'cavity'; k = (7+1/10)*sqrt(2)*pi; h = 1/6; degree = 5;
run(benchmark,degree,PREC,tol,iMax,iOut);
% BENCH WAVEGUIDE
iMax = 900; iOut = 30;
benchmark = 'waveguide'; k = 6*pi; h = 1/3.5; degree = 5;
run(benchmark,degree,PREC,tol,iMax,iOut);

% % BENCH FREE SPACE
% iMax = 1000; iOut = 50;
% benchmark = 'open'; k = 15*pi; h = 1/20; degree = 3;
% run(benchmark,degree,PREC,tol,iMax,iOut);
% % BENCH CAVITY
% iMax = 2000; iOut = 100;
% benchmark = 'cavity'; k = (7+1/10)*sqrt(2)*pi; h = 1/12; degree = 3;
% run(benchmark,degree,PREC,tol,iMax,iOut);
% % BENCH WAVEGUIDE
% iMax = 4000; iOut = 200;
% benchmark = 'waveguide'; k = 6*pi; h = 1/10; degree = 3;
% run(benchmark,degree,PREC,tol,iMax,iOut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,PREC,tol,iMax,iOut)
global k h

% Build mesh and dofManager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method CG (' benchmark ')']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);
[normErr] = computeNormError2D_CG(mesh, dofm, solA);

solP = computeSolProjL2_2D_CG(mesh, dofm);
normProjErr = computeNormError2D_CG(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(normErr, '%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(normProjErr, '%1.2e')]);
disp('---------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver CGNR']);
[resRedVec, resPhyVec, error1] = solverCGNRredu_CG(mesh, dofm, sysA, tol, iMax, iOut, @computeNormError2D_CG);

iterVec = (0:iOut:iMax)';
errorRef = normErr*ones(size(error1));

rezu1 = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2 = [iterVec resRedVec, resPhyVec, error1, errorRef];
name = sprintf('output/historyCGNR_CGred_%s_p%i_k%g_h%g.csv', benchmark, degree, k, h);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver GMRES']);
[resRedVec, resPhyVec, error2] = solverGMRESredu_CG(mesh, dofm, sysA, tol, iMax, iOut, @computeNormError2D_CG);

iterVec = (0:iOut:iMax)';
errorRef = normErr*ones(size(error2));

rezu1 = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2 = [iterVec resRedVec, resPhyVec, error2, errorRef];
name = sprintf('output/historyGMRES_CGred_%s_p%i_k%g_h%g.csv', benchmark, degree, k, h);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% errorNum = normErr*ones(size(error2));
% errorProj = normProjErr*ones(size(error2));
% 
% figure;
% hold off
% semilogy(iterVec,error1,'DisplayName','CGNR');
% hold on
% semilogy(iterVec,error2,'DisplayName','GMRES');
% semilogy(iterVec,errorNum,'k--','DisplayName','Numerical error');
% semilogy(iterVec,errorProj,'k:','DisplayName','Projection error');
% box on;
% grid on;
% title(['CG with condensation ' benchmark ' — k=' num2str(k/pi) 'pi — degree=' num2str(degree) ' — h=' num2str(h)]);
% legend('Location','southwest');
% xlabel('Iteration');
% ylabel('Relative error');
% axis([0 iMax 0.0025 1]);

end