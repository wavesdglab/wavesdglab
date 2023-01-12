%close all;
clear all;

headers2D;
global k;

degree = 3;
tau = 1;
theta = 1;
tol = 1e-100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE
% iMax = 1000; iOut = 50;
% benchmark = 'open'; k = 15*pi; h = 1/16;
% run(benchmark,degree,h,tau,theta,tol,iMax,iOut);
% benchmark = 'open'; k = 30*pi; h = 1/34;
% run(benchmark,degree,h,tau,theta,tol,iMax,iOut);

% BENCH CAVITY
iMax = 2000; iOut = 100;
benchmark = 'cavity'; k = (7+1/10)*sqrt(2)*pi; h = 1/10;
run(benchmark,degree,h,tau,theta,tol,iMax,iOut);
% benchmark = 'cavity'; k = (7+1/100)*sqrt(2)*pi; h = 1/15;
% run(benchmark,degree,h,tau,theta,tol,iMax,iOut);

% BENCH WAVEGUIDE
% iMax = 4000; iOut = 200;
% benchmark = 'waveguide'; k = 6*pi; h = 1/8;
% run(benchmark,degree,h,tau,theta,tol,iMax,iOut);
% benchmark = 'waveguide'; k = 12*pi; h = 1/17;
% run(benchmark,degree,h,tau,theta,tol,iMax,iOut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,h,tau,theta,tol,iMax,iOut)
global k;

% Build mesh and dofManager
mesh = benchmark2D(benchmark,h);
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method DG (' benchmark ')']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['    tau                 ' num2str(tau)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_DG(mesh, dofm, tau, theta);
[normErr, normErrU, normErrV, normSol, normSolU, normSolV] = computeNormError2D_DG(mesh, dofm, solA);

% [solP] = computeSolProjL2_2D_DG(mesh, dofm);
% [normProjErr, normProjErrU, normProjErrV] = computeNormError2D_DG(mesh, dofm, solP);
% 
% disp(['    L2-Norm Sol       ' num2str(normSol, '%1.2e') '  ' num2str(normSolU, '%1.2e') '  ' num2str(normSolV, '%1.2e')]);
% disp(['    L2-Norm ErrorSol  ' num2str(normErr, '%1.2e') '  ' num2str(normErrU, '%1.2e') '  ' num2str(normErrV, '%1.2e')]);
% disp(['    L2-Norm ErrorProj ' num2str(normProjErr, '%1.2e') '  ' num2str(normProjErrU, '%1.2e') '  ' num2str(normProjErrV, '%1.2e')]);
% disp('---------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% writeField_DG(dofm, mesh, solP, "output/mySol.pos", "mySol");
% system('gmsh output/mySol.pos');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver CGN']);
[resVec, error] = solverCGN_DG(mesh, dofm, sysA, tol, iMax, iOut);

iterVec = (0:iOut:iMax)';
errorRef = normErr*ones(size(error));

rezu1 = ["iter" "resVec" "error" "errorRef"];
rezu2 = [iterVec resVec, error, errorRef];
name = sprintf('output/historyCGN_DG_%s_p%i_k%g_h%g_tau%g+%gi.csv', benchmark, degree, k, h, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver GMRES']);
[resVec, error] = solverGMRES_DG(mesh, dofm, sysA, tol, iMax, iOut);

iterVec = (0:iOut:iMax)';
errorRef = normErr*ones(size(error));

rezu1 = ["iter" "resVec" "error" "errorRef"];
rezu2 = [iterVec resVec, error, errorRef];
name = sprintf('output/historyGMRES_DG_%s_p%i_k%g_h%g_tau%g+%gi.csv', benchmark, degree, k, h, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% hold off
% semilogy(iterVec,resVec  ,'-o','DisplayName','Relative residual');
% hold on
% semilogy(iterVec,error   ,'-o','DisplayName','Relative L2-error');
% semilogy(iterVec,errorRef,'k--','DisplayName','Relative L2-error (Ref)');
% box on;
% grid on;
% legend('Location','southwest');
% xlabel('Iteration');
% ylabel('Value');
% title(['DG ' benchmark ' — k=' num2str(k/pi) 'pi — h=' num2str(degree) ' — h=' num2str(h)]);
% axis([0 1000 1e-10 1]);

end