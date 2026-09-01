clear all;
close all;

benchmark = 'scattering_disk_penetrable';

% Parameters
global omega cAir cObj rhoAir rhoObj h
omega = 2.1*pi;
cAir = 1;
cObj = 0.5;
rhoAir = 1.33;
rhoObj = 0.6;
h = 1/10;

global Options
Options.Basis = 'Jacobi'; % Jacobi, Lobbato, Bernstein, Lagrange
degree = 3;
PREC = 0;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CG - Benchmark "' benchmark '"']);
disp(['---------------------------------------------------------']);
disp(['    omega               ' num2str(omega)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

% Compute numerical solution/error
[solA, sysA] = computeSolNum2D_CG_heterogeneous(mesh, dofm, PREC);
errorL2 = computeNormError2D_CG(mesh, dofm, solA);
disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);

% Compute projection solution/error
solP = computeSolProjL2_2D_CG(mesh, dofm);
errorProjL2 = computeNormError2D_CG(mesh, dofm, solP);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
global PML_HIDE; PML_HIDE = 1;
writeField2D(dofm, mesh, solA-solP, 'output/errNum.pos', "errNum");
system('gmsh output/mesh.msh output/solRef.pos output/solNum.pos output/errNum.pos output/mesh.msh&');
