%close all;
clear all;

global h
global omega eta1 eta2 k1 k2 c1 c2 rho1 rho2
global rho c eta k

degree = 3;
BASIS = 1;
PREC = 1;
A = 1;              % order of numerical fluxes
B = 2;              % order of transmission variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

benchmark = 'cavity_heterogeneous';
degree = 3; h = 1/10; omega = 15*pi;
rho1 = 1; c1 = 1; rho2 = 1; c2 = 1;
eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
run(benchmark,degree,BASIS,PREC,A,B);

benchmark = 'cavity_heterogeneous';
degree = 3; h = 1/10; omega = 15*pi;
rho1 = 1; c1 = 2; rho2 = 1; c2 = 0.8;
eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
run(benchmark,degree,BASIS,PREC,A,B);

benchmark = 'cavity_heterogeneous';
degree = 3; h = 1/10; omega = 15*pi;
rho1 = 1; c1 = 4; rho2 = 1; c2 = 1.6;
eta1 = rho1 * c1; eta2 = rho2 * c2; k1 = omega / c1; k2 = omega / c2;
run(benchmark,degree,BASIS,PREC,A,B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,BASIS,PREC,A,B)
global h
global omega eta1 eta2 k1 k2 c1 c2 rho1 rho2
global rho c eta k

mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);
setParameters(mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method CHDG (' benchmark ')']);
disp(['---------------------------------------------------------']);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_CHDG_ALL(mesh, dofm, PREC, A, B);
[normErr, ~, ~, normSol] = computeNormError2D_DG_ALL(mesh, dofm, solA);

solP = computeSolProjL2_2D_DG(mesh, dofm);
normProjErr = computeNormError2D_DG_ALL(mesh, dofm, solP);

disp(['    L2-Norm Sol       ' num2str(normSol, '%1.2e')]);
disp(['    L2-Norm ErrorSol  ' num2str(normErr, '%1.2e')]);
disp(['    L2-Norm ErrorProj ' num2str(normProjErr, '%1.2e')]);
disp('---------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sizeS = size(sysA.matS,1);
% disp(['    Size(S)             ' num2str(sizeS)]);
% nnzS = nnz(sysA.matS);
% disp(['    nnz(S)              ' num2str(nnzS)]);
% condS = condest(sysA.matS);
% disp(['    Condest(S)          ' num2str(condS, '%1.2e')]);
% condLocMin = min(condLoc);
% disp(['    CondMin(Loc)        ' num2str(condLocMin, '%1.2e')]);
% condLocMax = max(condLoc);
% disp(['    CondMax(Loc)        ' num2str(condLocMax, '%1.2e')]);
% disp(['---------------------------------------------------------']);
% 
% rezu1 = ["degree" "k" "h" "real(tau)" "imag(tau)" "sizeS" "nnzS" "normErr" "normProjErr" "normSol" "condS" "condLocMin" "condLocMax"];
% rezu2 = [degree, k, h, real(tau), imag(tau), sizeS, nnzS, normErr, normProjErr, normSol, condS, condLocMin, condLocMax];
% name = sprintf('output/statsCHDG_%s_p%i_k%g_h%g.csv', benchmark, degree, k, h);
% writematrix([rezu1 ; rezu2],name,'Delimiter','semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % % alpha = 1;
% % % % % matI = sparse(1:size(sysA.matS,1),1:size(sysA.matS,2),1);
% % % % % matIter = (1-alpha)*matI + alpha*(matI-sysA.matS);
% % % % % 
% % % % % % precL = inv(sysA.matGG);
% % % % % % precR = sparse(1:size(sysA.matGG,1),1:size(sysA.matGG,1),1);
% % % % % % precL = eigenvalGG\eigenvecGG';
% % % % % % precR = eigenvecGG;
% % % % % % matIter = precL*sysA.matS*precR;
% % % % % 
% % % % % % [eigenvecIter,eigenvalIter] = eigs(matIter,size(sysA.matS,1));
% % % % % % eigenvalIter                = diag(eigenvalIter);
% % % % % 
% % % % % eigenvalIter = eigs(matIter,size(sysA.matS,1));
% % % % % 
% % % % % rezu1 = ["real", "imag"];
% % % % % rezu2 = [real(eigenvalIter), imag(eigenvalIter)];
% % % % % name = sprintf('output/spectrumIter_%s_p%g_h%g.csv', benchmark, degree, h);
% % % % % writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');
% % % % % 
% % % % % disp(['    Min e.v. (Iter)     ' num2str(min(abs(eigenvalIter)))]);
% % % % % disp(['    Max e.v. (Iter)     ' num2str(max(abs(eigenvalIter)))]);
% % % % % % disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecIter))]);
% % % % % % disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecIter))]);
% % % % % disp('---------------------------------------------------------');

mat = sysA.matPinv*sysA.matS;

[~, eigenvalIter] = eigs(mat,size(mat,1));
eigenvalIter = 1 - diag(eigenvalIter);

rezu1 = ["real", "imag"];
rezu2 = [real(eigenvalIter), imag(eigenvalIter)];
name = sprintf('output/spectrumIter_%s_p%g_h%g.csv', benchmark, degree, h);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

disp(['    Min e.v. (Iter)     ' num2str(min(abs(eigenvalIter)))]);
disp(['    Max e.v. (Iter)     ' num2str(max(abs(eigenvalIter)))]);
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecIter))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecIter))]);
disp('---------------------------------------------------------');

benchmark = extractBefore(benchmark,'_');

figure();
hold off
scatter(real(eigenvalIter),imag(eigenvalIter),'DisplayName','Eigenvalues');
hold on
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k');
%plot(fovals(sysA.matS,100),'-b','DisplayName','Numerical range');
grid on; box on; axis equal;
title(['Benchmark "' benchmark '" — p=' num2str(degree) ' — h=' num2str(h)]);

end