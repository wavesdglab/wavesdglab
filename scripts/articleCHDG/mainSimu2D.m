clear;
close;

global omega k c1 c2 rho1 rho2 h1 h2

% Setup benchmark and parameters
benchmark = 'open_heterogeneous';
switch benchmark
    case 'open_heterogeneous'
        omega = 15*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 2; h1 = 1/16; h2 = 1/34;
        % omega = 15*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 1; h1 = 1/16; h2 = 1/34;
        tol = 1e-10; maxit = 1000; itout = 50;
    case 'disk_heterogeneous'
        omega = 10*pi; c1 = 1; c2 = 2/3; rho1 = 1; rho2 = 3/2; h1 = 0.06; h2 = 0.05;
        % omega = 10*pi; c1 = 1; c2 = 2/3; rho1 = 1; rho2 = 1; h1 = 0.1; h2 = 0.075;
        tol = 1e-10; maxit = 1000; itout = 100;
end

degree = 3;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);
Dlambda = 2*pi/k(1) * (sqrt(dofm.numDofTRI) - 1);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CHDG']);
disp(['---------------------------------------------------------']);
disp(['    h1                  ' num2str(h1)]);
disp(['    h2                  ' num2str(h2)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

% Compute numerical solution/error

% [solA, sysA] = computeSolNum2D_CHDG_sym(mesh, dofm, 1);
% errorL2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
% disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
% 
% [solA, sysA] = computeSolNum2D_CHDG_sym(mesh, dofm, 2);
% errorL2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
% disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);

[solA, sysA] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm);
errorL2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.8e')]);

[solA, sysA] = computeSolNum2D_CHDG_upw(mesh, dofm);
errorL2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.8e')]);

% [solA, sysA] = computeSolNum2D_HDG_sym(mesh, dofm, 0);
% errorL2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
% disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);

[solA, sysA] = computeSolNum2D_HDG_heterogeneous(mesh, dofm);
errorL2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.8e')]);

[solA, sysA] = computeSolNum2D_HDG_upw(mesh, dofm);
errorL2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.8e')]);

% Compute projection solution/error
solP = computeSolProjL2_2D_DG(mesh, dofm);
errorProjL2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solP);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
disp('---------------------------------------------------------');

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

% writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
% writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
% writeField2D(dofm, mesh, solA(1:mesh.numTri*3*dofm.numDofPerTRI)-solP, 'output/errNum.pos', "errNum");
% system('gmsh output/mesh.msh output/solRef.pos output/solNum.pos output/errNum.pos&');

% -------------------------------------------------------------------------
% Compute eigenvalues/eigenvectors
% -------------------------------------------------------------------------

% mat = sysA.matPinv*sysA.matS;
% [~, eigenval] = eigs(mat,size(mat,1));
% 
% eigenval = 1 - diag(eigenval);
% 
% 1 - max(abs(eigenval))
% 
% fprintf('Spectral radius = %.16f\n', max(abs(eigenval)));