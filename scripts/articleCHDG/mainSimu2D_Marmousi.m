clear;
close;

benchmark = 'geophysics_marmousi';

% Parameters
global omega nLambda
freq = 30;
omega = 2*pi*freq;
degree = 3;
nLambda = 10/(degree+1);

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

% Print coefficients
global c rho
writeCoef2D(mesh, c, 'output/velocity.pos', "Velocity [m/s]");
writeCoef2D(mesh, rho, 'output/density.pos', "Density");
system('output/mesh.msh gmsh output/velocity.pos output/density.pos&');

% -------------------------------------------------------------------------
% Compute solution
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CHDG - Benchmark "' benchmark '"']);
disp(['---------------------------------------------------------']);
disp(['    omega               ' num2str(omega/pi) '*pi']);
disp(['    nLambda             ' num2str(nLambda)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

disp('CHDG')
tic
[solCHDG, sysCHDG] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, 'SYM', 1);
toc
tic
writeField2D(dofm, mesh, solCHDG, 'output/solCHDG.pos', "CHDG");
toc

disp('HDG')
tic
[solHDG, sysHDG] = computeSolNum2D_HDG_heterogeneous(mesh, dofm, 'SYM');
toc
tic
writeField2D(dofm, mesh, solHDG, 'output/solHDG.pos', "HDG");
toc

disp('Difference')
tic
writeField2D(dofm, mesh, solCHDG-solHDG, 'output/solDiff.pos', "diff");
toc

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

system('gmsh output/mesh.msh output/velocity.pos output/density.pos output/solCHDG.pos output/solHDG.pos output/solDiff.pos&');
