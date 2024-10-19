clear all;
%close all;

benchmark = 'geophysics_marmousi';

% Parameters
global omega nLambda
freq = 5; %30
omega = 2*pi*freq;
degree = 3;
nLambda = 10/(degree+1);

%LREF = 1000;
%CREF = 1000;
%RHOREF = 1000;


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
%[solA, sysA] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm); % CHDG Sym-0
[solA, sysA] = computeSolNum2D_CHDG_upw(mesh, dofm); % CHDG Upw
disp('HDG')
%[solB, sysB] = computeSolNum2D_HDG_heterogeneous(mesh, dofm);  % HDG Sym-0
[solB, sysB] = computeSolNum2D_HDG_upw(mesh, dofm);  % HDG Upw

max(max(abs(solA-solB)))

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------
writeField2D(dofm, mesh, solA, 'output/solNumA.pos', "CHDG");
writeField2D(dofm, mesh, solB, 'output/solNumB.pos', "HDG");
system('gmsh output/mesh.msh output/velocity.pos output/density.pos output/solNumA.pos output/solNumB.pos&');
%writeField2D(dofm, mesh, solA-solB, 'output/diff.pos', "diff");
%system('gmsh output/mesh.msh output/solNumA.pos output/solNumB.pos output/diff.pos&');