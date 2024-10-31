clear;
close;

global omega c1 c2 rho1 rho2 h1 h2

degree = 3;
tol = 1e-100;
iMax = 1000;
iOut = 50;
PREC = 0;

% Setup benchmark and parameters

% benchmark = 'open_heterogeneous';
% omega = 15*pi; c1 = 1; c2 = 1;   rho1 = 1; rho2 = 1; h1 = 1/16; h2 = h1;
% omega = 30*pi; c1 = 1; c2 = 1;   rho1 = 1; rho2 = 1; h1 = 1/34; h2 = h1;
% omega = 15*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 2; h1 = 1/16; h2 = 1/34;
% omega = 15*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 1; h1 = 1/16; h2 = 1/34;

% benchmark = 'open_heterogeneous';
% omega = 3*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 1; h1 = 1/2; h2 = 1/4;
% omega = 3*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 1; h1 = 1/4; h2 = 1/8;
% omega = 3*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 1; h1 = 1/8; h2 = 1/16;
% omega = 3*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 1; h1 = 1/16; h2 = 1/32;
% 5.56e-02 8.61e-02 3.25e-02
% 4.28e-03 3.99e-02 2.82e-03
% 2.53e-04 2.79e-02 1.72e-04
% 1.71e-05 1.92e-02 1.16e-05

benchmark = 'disk_heterogeneous';
% omega = 16.5;  c1 = 1; c2 = 1;   rho1 = 1; rho2 = 1;   h1 = 0.04; h2 = h1;
% omega = 17;    c1 = 1; c2 = 1;   rho1 = 1; rho2 = 1;   h1 = 0.025; h2 = h1;
% omega = 36;    c1 = 1; c2 = 1;   rho1 = 1; rho2 = 1;   h1 = 0.065; h2 = h1;
% omega = 36.14; c1 = 1; c2 = 1;   rho1 = 1; rho2 = 1;   h1 = 0.055; h2 = h1;
% omega = 10*pi; c1 = 1; c2 = 2/3; rho1 = 1; rho2 = 3/2; h1 = 0.06;  h2 = 0.05;
% omega = 10*pi; c1 = 1; c2 = 2/3; rho1 = 1; rho2 = 1;   h1 = 0.1;   h2 = 0.075;

omega = 10*pi; c1 = 1; c2 = 2/3; rho1 = 1; rho2 = 1;   h1 = 0.1;   h2 = 0.1;
omega = 10*pi; c1 = 1; c2 = 2/3; rho1 = 1; rho2 = 1;   h1 = 0.1;   h2 = 0.05;
omega = 10*pi; c1 = 1; c2 = 2/3; rho1 = 1; rho2 = 1;   h1 = 0.1;   h2 = 0.025;

% Resonances modes of circular cavity
% J0 = @(x) besselj(0,x);
% for n = 1:10
%     z(n) = 2*fzero(J0,[(n-1) n]*pi);
% end
% 4.8097
% 11.0402
% 17.3075
% 23.5831
% 29.8618
% 36.1421

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

global k1 k2 eta1 eta2
disp(['---------------------------------------------------------']);
disp(['Benchmark ' benchmark '']);
disp(['---------------------------------------------------------']);
disp(['    k1                  ' num2str(k1)]);
disp(['    k2                  ' num2str(k2)]);
disp(['    eta1                ' num2str(eta1)]);
disp(['    eta2                ' num2str(eta2)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

% Compute numerical solution/error

% [solSym, sysSym] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, 'SYM', PREC);
% normErrSym = computeNormError2D_DG_heterogeneous(mesh, dofm, solSym);
% disp(['    L2-Error (CHDG SYM)   ' num2str(normErrSym,'%1.2e')]);

% [solSym2, sysSym2] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, 'SYM2', PREC);
% normErrSym2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solSym2);
% disp(['    L2-Error (CHDG SYM2)  ' num2str(normErrSym2,'%1.2e')]);

[solUpw, sysUpw] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, 'UPW', PREC);
normErrUpw = computeNormError2D_DG_heterogeneous(mesh, dofm, solUpw);
disp(['    L2-Error (CHDG UPW)   ' num2str(normErrUpw,'%1.2e')]);

[solUpw2, sysUpw2] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, 'UPW2', PREC);
normErrUpw2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solUpw2);
disp(['    L2-Error (CHDG UPW2)  ' num2str(normErrUpw2,'%1.2e')]);

% [solHdgSym, sysHdgSym] = computeSolNum2D_HDG_heterogeneous(mesh, dofm, 'SYM');
% normErrHdgSym = computeNormError2D_DG_heterogeneous(mesh, dofm, solHdgSym);
% disp(['    L2-Error (HDG SYM)    ' num2str(normErrHdgSym,'%1.2e')]);
% 
% [solHdgUpw, sysHdgUpw] = computeSolNum2D_HDG_heterogeneous(mesh, dofm, 'UPW');
% normErrHdgUpw = computeNormError2D_DG_heterogeneous(mesh, dofm, solHdgUpw);
% disp(['    L2-Error (HDG UPW)    ' num2str(normErrHdgUpw,'%1.2e')]);

solProj = computeSolProjL2_2D_DG(mesh, dofm);
normErrProj = computeNormError2D_DG_heterogeneous(mesh, dofm, solProj);
disp(['    L2-Error (projSol)  ' num2str(normErrProj,'%1.2e')]);
disp('---------------------------------------------------------');

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

% writeField2D(dofm, mesh, solSym, 'output/solNum.pos', "solNum");
% writeField2D(dofm, mesh, solProj, 'output/solRef.pos', "solRef");
% writeField2D(dofm, mesh, solSym(1:mesh.numTri*3*dofm.numDofPerTRI)-solProj, 'output/errNum.pos', "errNum");
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