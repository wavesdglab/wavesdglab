clear all;
%close all;

global omega c rho v0 h;

global Options
Options.Basis = 'Jacobi'; % Jacobi, Lobbato, Bernstein, Lagrange
degree = 3;

% Setup benchmark and parameters
benchmark = 'aeroacou_hydrointerf';
switch benchmark
    case 'open_convected'
        omega = 25*pi;
        h = 1/22;
        rho = 1;
        c = 1.5;
        k = omega / c;
        M = 1/6; % subsonic flow: 0<=M<1
        theta = pi/4;
        phi = 5*pi/4;
        v0 = [M*c*cos(theta), M*c*sin(theta)];
    case 'disk_convected'
        global theta v0d
        omega = 40; rho = 1; c = 1; theta = 0; v0d = 0.25; h = 1/20; % error 9.62e-03
        %omega = 40; rho = 1; c = 1; theta = 0; v0d = 0.5; h = 1/20; % error
        %omega = 40; rho = 1; c = 1; theta = 0; v0d = 0.75; h = 1/40; % error 1.23e-02
        v0 = v0d * [cos(theta), sin(theta)];
    case 'waveguide_convected'
        rho = 1;
        c = 1;
        omega = 5*pi;
        h = 1/11;
        eta = rho * c;
        k = omega / c;
        M = 0.5;   % subsonic flow: 0<=M<1
        v0 = [M*c, 0];
    case 'aeroacou_hydrointerf'
        rho = 1;
        c = 1.5;
        v0 = [0.5, 0.];
        lambda = 0.15;
        omega = 2*pi*v0(1)/lambda;
        degree = 5;
        h = 1/20;
end

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
