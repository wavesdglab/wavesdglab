close all;
clear;

global Options
Options.Basis = 'Jacobi';
Options.Error = 'Energy';

global omega c1 c2 rho1 rho2 h1 h2
degree = 3;
PREC = 0;

benchmark = 'open_heterogeneous';
% omega = 15*pi; c1 = 1; c2 = 1;   rho1 = 1; rho2 = 1; h1 = 1/16; h2 = h1;
% run(benchmark,degree,PREC);
% omega = 30*pi; c1 = 1; c2 = 1;   rho1 = 1; rho2 = 1; h1 = 1/34; h2 = h1;
% run(benchmark,degree,PREC);
omega = 15*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 2; h1 = 1/16; h2 = 1/34;
run(benchmark,degree,PREC);
omega = 15*pi; c1 = 1; c2 = 1/2; rho1 = 1; rho2 = 1; h1 = 1/16; h2 = 1/34;
run(benchmark,degree,PREC);

benchmark = 'disk_heterogeneous';
% omega = 16.5;  c1 = 1; c2 = 1;   rho1 = 1; rho2 = 1;   h1 = 0.04; h2 = h1;
% run(benchmark,degree,PREC);
% omega = 17;    c1 = 1; c2 = 1;   rho1 = 1; rho2 = 1;   h1 = 0.025; h2 = h1;
% run(benchmark,degree,PREC);
omega = 10*pi; c1 = 1; c2 = 2/3; rho1 = 1; rho2 = 3/2; h1 = 1/12;  h2 = 1/18;
run(benchmark,degree,PREC);
omega = 10*pi; c1 = 1; c2 = 2/3; rho1 = 1; rho2 = 1;   h1 = 1/12;  h2 = 1/18;
run(benchmark,degree,PREC);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,PREC)

global k1 k2 eta1 eta2

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

[solA, sysA] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, 'UPW', PREC);
normErrUpw = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
disp(['    L2-Error (CHDG UPW)   ' num2str(normErrUpw,'%1.2e')]);
%specRad = computeEigenval(sysA);
%name = sprintf('output/specRad_CHDG_UPW_%s_p%i_prec%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, PREC, k1, k2, eta1, eta2);
%writematrix(specRad, name, 'Delimiter', 'semi');

[solA, sysA] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, 'SYM', PREC);
normErrSym = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
disp(['    L2-Error (CHDG SYM)   ' num2str(normErrSym,'%1.2e')]);
%specRad = computeEigenval(sysA);
%name = sprintf('output/specRad_CHDG_SYM_%s_p%i_prec%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, PREC, k1, k2, eta1, eta2);
%writematrix(specRad, name, 'Delimiter', 'semi');

[solA, sysA] = computeSolNum2D_CHDG_heterogeneous(mesh, dofm, 'SYM2', PREC);
normErrSym2 = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
disp(['    L2-Error (CHDG SYM2)  ' num2str(normErrSym2,'%1.2e')]);
%specRad = computeEigenval(sysA);
%name = sprintf('output/specRad_CHDG_SYM2_%s_p%i_prec%i_k_%g_%g_eta_%g_%g.csv', benchmark, degree, PREC, k1, k2, eta1, eta2);
%writematrix(specRad, name, 'Delimiter', 'semi');

[solA, sysA] = computeSolNum2D_HDG_heterogeneous(mesh, dofm, 'SYM');
normErrHdgSym = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
disp(['    L2-Error (HDG SYM)    ' num2str(normErrHdgSym,'%1.2e')]);

[solA, sysA] = computeSolNum2D_HDG_heterogeneous(mesh, dofm, 'UPW');
normErrHdgUpw = computeNormError2D_DG_heterogeneous(mesh, dofm, solA);
disp(['    L2-Error (HDG UPW)    ' num2str(normErrHdgUpw,'%1.2e')]);

solProj = computeSolProjL2_2D_DG(mesh, dofm);
normErrProj = computeNormError2D_DG_heterogeneous(mesh, dofm, solProj);
disp(['    L2-Error (projSol)  ' num2str(normErrProj,'%1.2e')]);
disp('---------------------------------------------------------');

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

% writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
% writeField2D(dofm, mesh, solProj, 'output/solRef.pos', "solRef");
% writeField2D(dofm, mesh, solA(1:mesh.numTri*3*dofm.numDofPerTRI)-solProj, 'output/errNum.pos', "errNum");
% system('gmsh output/mesh.msh output/solRef.pos output/solNum.pos output/errNum.pos&');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Compute eigenvalues/eigenvectors
% -------------------------------------------------------------------------

function specRad = computeEigenval(sysA)

tic
[~, eigenval] = eigs(sysA.matPinv*sysA.matS,size(sysA.matS,1));
eigenval = 1 - diag(eigenval);
specRad = 1-max(abs(eigenval));
fprintf('Spectral radius = %.16f\n', 1-max(abs(eigenval)));
toc

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%