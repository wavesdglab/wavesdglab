close all;
clear all;

global k h R_disk L l L_PML;

% Setup benchmark and parameters
benchmark = 'scattering';
k = 10;
h = 0.05;
L = 1.1;
l = L;
R_disk = 1;
L_PML = 0.3;
degree = 1;
PREC = 0;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

disp(['---------------------------------------------------------']);
disp(['Method CG - Benchmark "' benchmark '"']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

%[solA, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);
[solA, sysA] = computeSolNum2D_PML_CG(mesh, dofm, PREC);
size(solA)
errorL2 = computeNormError2D_PML_CG(mesh, dofm, solA);
disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);

solP = computeSolProjL2_2D_CG(mesh, dofm);
errorL2 = computeNormError2D_CG(mesh, dofm, solP);
disp(['    L2-Error (refSol)   ' num2str(errorL2,'%1.2e')]);

writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
writeField2D(dofm, mesh, solA-solP, 'output/errNum.pos', "errNum");
system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

