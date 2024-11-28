clear all;
%close all;

global omega k eta c rho h v0 M l;
   
% Setup benchmark and parameters
benchmark = 'waveguide_convected';
switch benchmark
        case 'waveguide_convected'
        omega = 5*pi;
        h = 1/25;
        l = 1;
        tol = 1e-10; maxit = 1000; itout = 50;
        rho = 1;
        c = 1;  
        eta = rho * c;
        k = omega / c;
        M = 0.05;           % subsonic flow: 0<=M<1
        v0 = [M*c, 0];
end
degree = 3;
BASIS = 0;
PREC = 1;           

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);
Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

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

% Compute numerical and projection solutions
disp('---------------------------------------------------------');
disp('    HDG   ');
[solA, sysA] = computeSolNum2D_HDG_convected(mesh, dofm, PREC);
disp('---------------------------------------------------------');
disp('    CHDG   ');
[solB, sysB] = computeSolNum2D_CHDG_convected(mesh, dofm, PREC);
disp('---------------------------------------------------------');
disp('    Projection   ');
solP = computeSolProjL2_2D_DG(mesh, dofm);
disp('---------------------------------------------------------');

% Compute numerical and projection errors
errorL2_A = computeNormError2D_DG_convected(mesh, dofm, solA);
errorL2_B = computeNormError2D_DG_convected(mesh, dofm, solB);
errorProjL2 = computeNormError2D_DG_convected(mesh, dofm, solP);

% Display numerical and projection errors
disp('---------------------------------------------------------');
disp(['    L2-Error (numSol)  HDG          ' num2str(errorL2_A,'%1.2e')]);
disp(['    L2-Error (numSol)  CHDG         ' num2str(errorL2_B,'%1.2e')]);
disp(['    L2-Error (projSol) Projection   ' num2str(errorProjL2,'%1.2e')]);
disp('---------------------------------------------------------');

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
writeField2D(dofm, mesh, solA(1:mesh.numTri*3*dofm.numDofPerTRI)-solP, 'output/errNum.pos', "errNum");
system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

writeField2D(dofm, mesh, solB, 'output/solNum2.pos', "solNum");
writeField2D(dofm, mesh, solP, 'output/solRef2.pos', "solRef");
writeField2D(dofm, mesh, solB(1:mesh.numTri*3*dofm.numDofPerTRI)-solP, 'output/errNum2.pos', "errNum");
system('gmsh output/solRef2.pos output/solNum2.pos output/errNum2.pos&');