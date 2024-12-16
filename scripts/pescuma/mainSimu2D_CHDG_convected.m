clear all;
%close all;

global omega k eta c rho h v0 M theta phi;

% H=[1/10 1/15 1/20 1/25];
% Degree=[1 2 3 4];

% for i=1:size(H,2)
%     h=H(i);
%     degree = Degree(i);

% Setup benchmark and parameters
benchmark = 'open_convected';
switch benchmark
    case 'waveguide_convected'
        omega = 5*pi;
        h = 1/10;
        tol = 1e-10; maxit = 1000; itout = 50;
        rho = 1;
        c = 1;  
        eta = rho * c;
        k = omega / c;
        M = 0.1;           % subsonic flow: 0<=M<1
        v0 = [M*c, 0];
    case 'open_convected'
        omega = sqrt(2)*10;
        h = 1/15;
        tol = 1e-10; maxit = 1000; itout = 50;
        rho = 1;
        c = 1;  
        eta = rho * c;
        k = omega / c;
        M = 0.3;           % subsonic flow: 0<=M<1
        theta = pi/6;
        phi = pi/3;
        v0 = [M*c*cos(theta), M*c*sin(theta)];
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
disp('    DG   ');
[solA, sysA] = computeSolNum2D_DG_convected(mesh, dofm, PREC);
disp('---------------------------------------------------------');

disp('---------------------------------------------------------');
disp('    HDG   ');
[solB, sysB] = computeSolNum2D_HDG_convected(mesh, dofm, PREC);
disp('---------------------------------------------------------');

% disp('---------------------------------------------------------');
% disp('    CHDG   ');
% [solC, sysC] = computeSolNum2D_CHDG_convected(mesh, dofm, PREC);
% disp('---------------------------------------------------------');

disp('---------------------------------------------------------');
disp('    Projection   ');
solP = computeSolProjL2_2D_DG(mesh, dofm);
disp('---------------------------------------------------------');

% Compute numerical and projection errors
errorL2_A = computeNormError2D_DG_convected(mesh, dofm, solA);
errorL2_B = computeNormError2D_DG_convected(mesh, dofm, solB);
% errorL2_C = computeNormError2D_DG_convected(mesh, dofm, solC);
errorProjL2 = computeNormError2D_DG_convected(mesh, dofm, solP);

% Display numerical and projection errors
disp('---------------------------------------------------------');
disp(['    L2-Error (numSol)  DG          ' num2str(errorL2_A,'%1.2e')]);
disp(['    L2-Error (numSol)  HDG          ' num2str(errorL2_B,'%1.2e')]);
% disp(['    L2-Error (numSol)  CHDG         ' num2str(errorL2_C,'%1.2e')]);
disp(['    L2-Error (projSol) Projection   ' num2str(errorProjL2,'%1.2e')]);
disp('---------------------------------------------------------');

% err(i)=errorL2_A;

% end

% figure(1)
% semilogy(H,err);
% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------
writeField2D(dofm, mesh, solP, 'output/solRef.pos', "Ref");

writeField2D(dofm, mesh, solA, 'output/solDG.pos', "DG");
writeField2D(dofm, mesh, solA(1:mesh.numTri*3*dofm.numDofPerTRI)-solP, 'output/errNumDG.pos', "errNumDG");
system('gmsh output/mesh.msh output/solDG.pos output/solRef.pos output/errNumDG.pos&');

writeField2D(dofm, mesh, solB, 'output/solHDG.pos', "HDG");
writeField2D(dofm, mesh, solB(1:mesh.numTri*3*dofm.numDofPerTRI)-solP, 'output/errNumHDG.pos', "errNumHDG");
system('gmsh output/mesh.msh output/solHDG.pos output/solRef.pos output/errNumHDG.pos&');

% writeField2D(dofm, mesh, solC, 'output/solCHDG.pos', "CHDG");
% writeField2D(dofm, mesh, solC(1:mesh.numTri*3*dofm.numDofPerTRI)-solP, 'output/errNumCHDG.pos', "errNumCHDG");
% system('gmsh output/mesh.msh output/solCHDG.pos output/solRef.pos output/errNumCHDG.pos&');

diff = solA - solB;
norm_infinity_diff = max(max(abs(solA-solB)))
writeField2D(dofm, mesh, diff, 'output/soldiff.pos', "diff");
system('gmsh output/mesh.msh output/soldiff.pos&');