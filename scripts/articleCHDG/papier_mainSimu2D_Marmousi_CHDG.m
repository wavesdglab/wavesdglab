clear all;
close all;

benchmark = 'geophysics_marmousi';

% Parameters
global omega nLambda
freq = 5; %30
omega = 2*pi*freq;
degree = 3;
nLambda = 10/(degree+1);
PREC = 0;
A = 1;
B = 1;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

% Print coefficients
global cArray rhoArray
writeCoef2D(mesh, cArray, 'output/velocity.pos', "Velocity [m/s]");
writeCoef2D(mesh, rhoArray, 'output/density.pos', "Density");
system('gmsh output/velocity.pos output/density.pos&');

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CHDG - Benchmark "' benchmark '"']);
disp(['---------------------------------------------------------']);
disp(['    omega               ' num2str(omega)]);
disp(['    nLambda             ' num2str(nLambda)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

% Compute numerical solution/error
% [solA, sysA] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, PREC, A, B);
% [solA, sysA] = computeSolNum2D_CHDG_heterogeneous_2(mesh, dofm, PREC); % CHDG with 0th-order symmetric fluxes (only blocks I and G)
[solB, sysB] = computeSolNum2D_HDG_heterogeneous_2(mesh, dofm, PREC);    % HDG with 0th-order symmetric fluxes (only blocks I and G)

% sysB.matGGinv*sysB.matGG

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

% writeField2D(dofm, mesh, solA, 'output/solNumA.pos', "solNum");
% system('gmsh output/solNumA.pos&');

writeField2D(dofm, mesh, solB, 'output/solNumB.pos', "solNum");
system('gmsh output/solNumB.pos&');