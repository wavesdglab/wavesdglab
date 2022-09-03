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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build mesh and dofManager
mesh = benchmark2D(benchmark,h);
mesh = buildMeshConnectivity(mesh);
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

[solA, sysA] = computeSolNum2D_CG(mesh, dofm);
[errorL2] = computeNormError2D_CG(mesh, dofm, solA);

[solP, sysP] = computeSolProjL2_2D_CG(mesh, dofm);
[errorProjL2] = computeNormError2D_CG(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(errorL2, '%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2, '%1.2e')]);
disp('---------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% writeField_CG(dofm, mesh, solP, "output/mySol.pos", "mySol");
% system('gmsh output/mySol.pos');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver CGN']);
tol = 1e-10;
iMax = 1000;
iOut = 50;
[resRedVec, resPhyVec, error] = solverCGNredu_CG(mesh, dofm, sysA, tol, iMax, iOut);

iterVec = (0:iOut:iMax)';
errorRef = errorL2*ones(size(error));

rezu1 = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2 = [iterVec resRedVec, resPhyVec, error, errorRef];
name = sprintf('output/historyCGN_CG_%s_P%i_k%g_h%g.csv', benchmark, degree, k, h);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['--- Solver GMRES']);
tol = 1e-10;
iMax = 1000;
iOut = 50;
[resRedVec, resPhyVec, error] = solverGMRESredu_CG(mesh, dofm, sysA, tol, iMax, iOut);

iterVec = (0:iOut:iMax)';
errorRef = errorL2*ones(size(error));

rezu1 = ["iter" "resRed" "resPhy" "error" "errorRef"];
rezu2 = [iterVec resRedVec, resPhyVec, error, errorRef];
name = sprintf('output/historyGMRES_CG_%s_P%i_k%g_h%g.csv', benchmark, degree, k, h);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% hold off
% semilogy(iterVec,resPhyVec   ,'-o','DisplayName','Relative residual (Phy)');
% hold on
% semilogy(iterVec,resRedVec   ,'-x','DisplayName','Relative residual (Red)');
% semilogy(iterVec,error       ,'-o','DisplayName','Relative L2-error');
% semilogy(iterVec,errorRef    ,'k--','DisplayName','Relative L2-error (Ref)');
% box on;
% grid on;
% legend('Location','southwest');
% xlabel('Iteration');
% ylabel('Value');
% %axis([0 1000 1e-7 1]);