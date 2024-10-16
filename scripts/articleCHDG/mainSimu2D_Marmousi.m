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

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

% Print coefficients
global c rho
writeCoef2D(mesh, c, 'output/velocity.pos', "Velocity [m/s]");
writeCoef2D(mesh, rho, 'output/density.pos', "Density");
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
disp('CHDG')
[solA, sysA] = computeSolNum2D_CHDG_Marmousi(mesh, dofm, PREC); % CHDG Sym-0
disp('HDG')
[solB, sysB] = computeSolNum2D_HDG_Marmousi(mesh, dofm, PREC);  % HDG Sym-0

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------
writeField2D(dofm, mesh, solA, 'output/solNumA.pos', "CHDG");
writeField2D(dofm, mesh, solB, 'output/solNumB.pos', "HDG");
writeField2D(dofm, mesh, solA-solB, 'output/diff.pos', "diff");
system('gmsh output/solNumA.pos output/solNumB.pos output/diff.pos&');