clear all;
%close all;

global omega c rho v0 h;

global Options
Options.Basis = 'Jacobi'; % Jacobi, Lobbato, Bernstein, Lagrange
degree = 3;

% Setup benchmark and parameters
benchmark = 'aeroacou_hydrointerf';
rho = 1;
c = 1.5;
v0 = [0.5, 0.];
lambda = 0.15;
omega = 2*pi*v0(1)/lambda;
degree = 5;
h = 1/20;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

% Compute numerical solution/error
solNum = computeSolNum2D_CHDG_convected(mesh, dofm, 1);
solP = solNum((1:dofm.numDofTRI)+0*dofm.numDofTRI);
solUx = solNum((1:dofm.numDofTRI)+1*dofm.numDofTRI);
solUy = solNum((1:dofm.numDofTRI)+2*dofm.numDofTRI);

writeField2D(dofm, mesh, solP, 'output/solP.pos', "solP");
writeField2D(dofm, mesh, solUx, 'output/solUx.pos', "solUx");
writeField2D(dofm, mesh, solUy, 'output/solUy.pos', "solUy");

system('gmsh output/mesh.msh output/solP.pos output/solUx.pos output/solUy.pos&');
