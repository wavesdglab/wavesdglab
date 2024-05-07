clear all;
close all;

% benchmark = 'geophysics_BPmodel';
benchmark = 'geophysics_marmousi';

% Parameters
global omega nLambda
omega = 100*pi;  % FIXME: THERE IS SOMETHING NOT CLEAR FOR THE CHOICE OF 'LC'
nLambda = 0.2;

degree = 3;
PREC = 0;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

% Print coefficients
global cArray rhoArray
writeCoef2D(mesh, cArray, 'output/velocity.pos', "Velocity [m/s]");
writeCoef2D(mesh, rhoArray, 'output/density.pos', "Density");
system('gmsh output/velocity.pos output/density.pos&');

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

writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
system('gmsh output/solNum.pos&');
