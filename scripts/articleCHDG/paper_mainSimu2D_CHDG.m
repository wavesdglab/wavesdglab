clear all;
%close all;

global omega k eta1 eta2 k1 k2 c1 c2 rho1 rho2 h1 h2
   
% Setup benchmark and parameters
benchmark = 'disk_heterogeneous';
switch benchmark
    case 'open_heterogeneous'
        omega = 15*pi; %15*pi;
        h1 = 1/16;
        h2 = 1/16;
        tol = 1e-10; maxit = 1000; itout = 50;
        rho1 = 1;
        c1 = 1;  
        rho2 = 1;
        c2 = 1; 
        eta1 = rho1 * c1;
        eta2 = rho2 * c2;
        k1 = omega / c1;
        k2 = omega / c2;
    case 'disk_heterogeneous'
        omega = 36; %10*pi;
        h1 = 0.065;
        h2 = 0.065;
        tol = 1e-10; maxit = 1000; itout = 100;
        rho1 = 1;
        c1 = 1; 
        rho2 = 1;
        c2 = 1;
        eta1 = rho1 * c1;
        eta2 = rho2 * c2;
        k1 = omega / c1;
        k2 = omega / c2;
end
degree = 3;
BASIS = 0;
PREC = 1;
A = 2;              % order of numerical fluxes: A=1 for 0th order, A=2 for 2nd order            

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

% Compute numerical solution
[solA, sysA] = computeSolNum2D_CHDG_sym(mesh, dofm, PREC, A);
% [solA, sysA] = computeSolNum2D_CHDG_upw(mesh, dofm, PREC);

% Compute numerical error
errorL2_A = computeNormError2D_DG_ALL(mesh, dofm, solA);
% errorL2_B = computeNormError2D_DG_ALL(mesh, dofm, solB)
% errorL2_C = computeNormError2D_DG_ALL(mesh, dofm, solC)
errorL2 = errorL2_A

% Compute projection solution
solP = computeSolProjL2_2D_DG(mesh, dofm);
% errorProjL2 = computeNormError2D_DG_ALL(mesh, dofm, solP);
% 
% disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
% disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
% disp('---------------------------------------------------------');

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
writeField2D(dofm, mesh, solA(1:mesh.numTri*3*dofm.numDofPerTRI)-solP, 'output/errNum.pos', "errNum");
system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

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