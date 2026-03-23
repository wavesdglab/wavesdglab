clear all;
%close all;

global omega c rho h v0 v0d;

global Options
Options.Basis = 'Jacobi'; % Jacobi, Lobbato, Bernstein, Lagrange
Options.Error = 'Energy';
degree = 3;
BASIS = 0;
PREC = 0;

% Setup benchmark and parameters
benchmark = 'disk_convected';
switch benchmark
    case 'open_convected'
        omega = 25*pi;
        h = 1/22;
        tol = 1e-10; maxit = 1000; itout = 50;
        rho = 1;
        c = 1.5;
        k = omega / c;
        M = 1/6; % subsonic flow: 0<=M<1
        theta = pi/4;
        phi = 5*pi/4;
        v0 = [M*c*cos(theta), M*c*sin(theta)];
    case 'disk_convected'
        %omega = 40; rho = 1; c = 1; theta = 0; v0d = 0.25; h = 1/20; % error 9.62e-03
        omega = 40; rho = 1; c = 1; theta = 0; v0d = 0.5; h = 1/20; % error   
        %omega = 40; rho = 1; c = 1; theta = 0; v0d = 0.75; h = 1/40; % error 1.23e-02
        v0 = v0d * [cos(theta), sin(theta)];
        tol = 1e-10; maxit = 1000; itout = 50;
        Dlambda = 2*pi*c/omega / h * (degree+1);
    case 'waveguide_convected'
        omega = 5*pi;
        h = 1/11;
        tol = 1e-10; maxit = 1000; itout = 50;
        rho = 1;
        c = 1;
        eta = rho * c;
        k = omega / c;
        M = 0.5;   % subsonic flow: 0<=M<1
        v0 = [M*c, 0];
end

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);
if(isempty(Dlambda))
    Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);
end

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CHDG']);
disp(['---------------------------------------------------------']);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

solP = computeSolProjL2_2D_DG(mesh, dofm);
[normProjErr, normSol] = computeNormError2D_DG_convected(mesh, dofm, solP);
writeField2D(dofm, mesh, solP, 'output/solProj.pos', "SolProj");

disp(['    L2-Norm Sol       ' num2str(normSol, '%1.2e')]);
disp(['    L2-Norm ErrorProj ' num2str(normProjErr, '%1.2e')]);

% [solA, sysA] = computeSolNum2D_DG_convected(mesh, dofm);
% normErrA = computeNormError2D_DG_convected(mesh, dofm, solA);
% disp(['    L2-Norm ErrorDG   ' num2str(normErrA, '%1.2e')]);

% [solB, sysB] = computeSolNum2D_HDG_convected(mesh, dofm, PREC);
% errorL2_B = computeNormError2D_DG_convected(mesh, dofm, solB);
% disp(['    L2-Norm ErrorHDG  ' num2str(normErrA, '%1.2e')]);

[solC, sysC] = computeSolNum2D_CHDG_convected(mesh, dofm, PREC);
normErrC = computeNormError2D_DG_convected(mesh, dofm, solC);
writeField2D(dofm, mesh, solC, 'output/solNum.pos', "solNum");
writeField2D(dofm, mesh, solC-solP, 'output/errNum.pos', "errNum");

disp(['    L2-Norm ErrorCHDG ' num2str(normErrC, '%1.2e')]);
disp('---------------------------------------------------------');

system('gmsh output/mesh.msh output/solProj.pos output/solNum.pos output/errNum.pos&');
