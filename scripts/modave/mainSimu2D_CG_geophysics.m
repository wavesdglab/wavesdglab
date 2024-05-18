clear all;
close all;

% benchmark = 'geophysics_BPmodel'; freq = 20;
benchmark = 'geophysics_marmousi'; freq = 50;

% Parameters
global omega nLambda
omega = 2*pi*freq;
degree = 3;
nLambda = 10/(degree+1);
PREC = 0;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

% Print coefficients
global cArray rhoArray
writeCoef2D(mesh, cArray, 'output/velocity.pos', "Velocity [m/s]");
writeCoef2D(mesh, rhoArray, 'output/density.pos', "Density");
system('gmsh output/mesh.msh output/velocity.pos output/density.pos&');

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CG - Benchmark "' benchmark '"']);
disp(['---------------------------------------------------------']);
disp(['    omega               ' num2str(omega)]);
disp(['    nLambda             ' num2str(nLambda)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

% Compute numerical solution/error
[solA, sysA] = computeSolNum2D_CG_heterogeneous(mesh, dofm, PREC);

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

writeField2D(dofm, mesh, solA, 'output/solution.pos', "Solution");
system('gmsh output/mesh.msh output/velocity.pos output/density.pos output/solution.pos output/mesh.msh&');