%close all;
clear all;

headers2D;
global k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE
% benchmark = 'open'; degree = 1; k = 10*pi; h = 1/32;
benchmark = 'open'; degree = 3; k = 10*pi; h = 1/8;

% BENCH CAVITY
% benchmark = 'cavity'; degree = 1; k = (3+1/8)*sqrt(2)*pi; h = 1/32;
% benchmark = 'cavity'; degree = 3; k = (5+1/8)*sqrt(2)*pi; h = 1/8;

% BENCH WAVEGUIDE
% benchmark = 'waveguide'; degree = 1; k = 2*pi; h = 1/16;
% benchmark = 'waveguide'; degree = 3; k = 6*pi; h = 1/8;

tau = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build mesh and dofManager
mesh = benchmark2D(benchmark,h);
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method HDG (' benchmark ')']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['    tau                 ' num2str(tau)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_HDG(mesh, dofm, tau);
[errorL2] = computeNormError2D_DG(mesh, dofm, solA);

[solP, sysP] = computeSolProjL2_2D_DG(mesh, dofm);
[errorProjL2] = computeNormError2D_DG(mesh, dofm, solP);

[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solA);
[errorPostL2] = computeNormError2D_DG(mesh, dofmPost, solApost);

[solPpost, sysPpost] = computeSolProjL2_2D_DG(mesh, dofmPost);
[errorProjPostL2] = computeNormError2D_DG(mesh, dofmPost, solPpost);

disp(['    L2-Error (numSol)   ' num2str(errorL2, '%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2, '%1.2e')]);
disp(['    L2-Error (numPost)  ' num2str(errorPostL2, '%1.2e')]);
disp(['    L2-Error (projPost) ' num2str(errorProjPostL2, '%1.2e')]);
disp('---------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% writeField_DG(dofm, mesh, solP, "output/mySol.pos", "mySol");
% system('gmsh output/mySol.pos');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(['--- Solver Richardson']);
% tol = 1e-10;
% iMax = 100000;
% iOut = 5000;
% alpha = 1;
% [resRedVec, resPhyVec, error, errorPost] = solverRichardson_DG(mesh, dofm, sysA, tol, iMax, iOut, alpha);
% 
% iterVec = (0:iOut:iMax)';
% errorRef = errorL2*ones(size(error));
% errorRefPost = errorPostL2*ones(size(error));
% 
% rezu1 = ["iter" "resRed" "resPhy" "error" "errorPost" "errorRef" "errorRefPost"];
% rezu2 = [iterVec resRedVec, resPhyVec, error, errorPost, errorRef, errorRefPost];
% name = sprintf('output/historyRich_HDG_%s_P%i_k%g_h%g_tau%g+%gi_alpha%g.csv', benchmark, degree, k, h, real(tau), imag(tau), alpha);
% writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver CGN']);
tol = 1e-10;
iMax = 1000;
iOut = 50;
[resRedVec, resPhyVec, error, errorPost] = solverCGNredu_DG(mesh, dofm, sysA, tol, iMax, iOut);

iterVec = (0:iOut:iMax)';
errorRef = errorL2*ones(size(error));
errorRefPost = errorPostL2*ones(size(error));

rezu1 = ["iter" "resRed" "resPhy" "error" "errorPost" "errorRef" "errorRefPost"];
rezu2 = [iterVec resRedVec, resPhyVec, error, errorPost, errorRef, errorRefPost];
name = sprintf('output/historyCGN_HDG_%s_P%i_k%g_h%g_tau%g+%gi.csv', benchmark, degree, k, h, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver GMRES']);
tol = 1e-10;
iMax = 1000;
iOut = 50;
[resRedVec, resPhyVec, error, errorPost] = solverGMRESredu_DG(mesh, dofm, sysA, tol, iMax, iOut);

iterVec = (0:iOut:iMax)';
errorRef = errorL2*ones(size(error));
errorRefPost = errorPostL2*ones(size(error));

rezu1 = ["iter" "resRed" "resPhy" "error" "errorPost" "errorRef" "errorRefPost"];
rezu2 = [iterVec resRedVec, resPhyVec, error, errorPost, errorRef, errorRefPost];
name = sprintf('output/historyGMRES_HDG_%s_P%i_k%g_h%g_tau%g+%gi.csv', benchmark, degree, k, h, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% hold off
% semilogy(iterVec,resPhyVec   ,'-o','DisplayName','Relative residual (Phy)');
% hold on
% semilogy(iterVec,resRedVec   ,'-x','DisplayName','Relative residual (Red)');
% semilogy(iterVec,error       ,'-o','DisplayName','Relative L2-error');
% semilogy(iterVec,errorPost   ,'-x','DisplayName','Relative L2-error with PostPro');
% semilogy(iterVec,errorRef    ,'k--','DisplayName','Relative L2-error (Ref)');
% semilogy(iterVec,errorRefPost,'k--','DisplayName','Relative L2-error with PostPro (Ref)');
% box on;
% grid on;
% legend('Location','southwest');
% xlabel('Iteration');
% ylabel('Value');
% %axis([0 1000 1e-4 1]);