%close all;
clear all;

global k;

tau = 1;
prec = 10;
degree = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE
benchmark = 'open'; k = 15*pi; h = 1/16;
run(benchmark,degree,h,tau,prec);
benchmark = 'open'; degree = 3; k = 15*pi; h = 1/34;
run(benchmark,degree,h,tau,prec);
benchmark = 'open'; degree = 3; k = 30*pi; h = 1/34;
run(benchmark,degree,h,tau,prec);

% BENCH CAVITY
benchmark = 'cavity'; k = (7+1/10)*sqrt(2)*pi; h = 1/10;
run(benchmark,degree,h,tau,prec);
benchmark = 'cavity'; degree = 3; k = (7+1/10)*sqrt(2)*pi; h = 1/15;
run(benchmark,degree,h,tau,prec);
benchmark = 'cavity'; degree = 3; k = (7+1/100)*sqrt(2)*pi; h = 1/15;
run(benchmark,degree,h,tau,prec);

% BENCH WAVEGUIDE
benchmark = 'waveguide'; k = 6*pi; h = 1/8;
run(benchmark,degree,h,tau,prec);
benchmark = 'waveguide'; degree = 3; k = 6*pi; h = 1/17;
run(benchmark,degree,h,tau,prec);
benchmark = 'waveguide'; degree = 3; k = 12*pi; h = 1/17;
run(benchmark,degree,h,tau,prec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,h,tau,prec)
global k;

mesh = setupBenchmark2D(benchmark,h);
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method CHDG (' benchmark ')']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['    tau                 ' num2str(tau)]);
disp(['---------------------------------------------------------']);

[solA, sysA, condLoc] = computeSolNum2D_CHDG(mesh, dofm, tau, prec);
[normErr, ~, ~, normSol] = computeNormError2D_DG(mesh, dofm, solA);

solP = computeSolProjL2_2D_DG(mesh, dofm);
normProjErr = computeNormError2D_DG(mesh, dofm, solP);

disp(['    L2-Norm Sol       ' num2str(normSol, '%1.2e')]);
disp(['    L2-Norm ErrorSol  ' num2str(normErr, '%1.2e')]);
disp(['    L2-Norm ErrorProj ' num2str(normProjErr, '%1.2e')]);
disp('---------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sizeS = size(sysA.matS,1);
disp(['    Size(S)             ' num2str(sizeS)]);
nnzS = nnz(sysA.matS);
disp(['    nnz(S)              ' num2str(nnzS)]);
condS = condest(sysA.matS);
disp(['    Condest(S)          ' num2str(condS, '%1.2e')]);
condLocMin = min(condLoc);
disp(['    CondMin(Loc)        ' num2str(condLocMin, '%1.2e')]);
condLocMax = max(condLoc);
disp(['    CondMax(Loc)        ' num2str(condLocMax, '%1.2e')]);
disp(['---------------------------------------------------------']);

rezu1 = ["degree" "k" "h" "real(tau)" "imag(tau)" "sizeS" "nnzS" "normErr" "normProjErr" "normSol" "condS" "condLocMin" "condLocMax"];
rezu2 = [degree, k, h, real(tau), imag(tau), sizeS, nnzS, normErr, normProjErr, normSol, condS, condLocMin, condLocMax];
name = sprintf('output/statsCHDG_%s_p%i_k%g_h%g.csv', benchmark, degree, k, h);
writematrix([rezu1 ; rezu2],name,'Delimiter','semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alpha = 1;
% matI = sparse(1:size(sysA.matS,1),1:size(sysA.matS,2),1);
% matIter = (1-alpha)*matI + alpha*(matI-sysA.matS);
% 
% % precL = inv(sysA.matGG);
% % precR = sparse(1:size(sysA.matGG,1),1:size(sysA.matGG,1),1);
% % precL = eigenvalGG\eigenvecGG';
% % precR = eigenvecGG;
% % matIter = precL*sysA.matS*precR;
% 
% % [eigenvecIter,eigenvalIter] = eigs(matIter,size(sysA.matS,1));
% % eigenvalIter                = diag(eigenvalIter);
% 
% eigenvalIter = eigs(matIter,size(sysA.matS,1));
% 
% rezu1 = ["real", "imag"];
% rezu2 = [real(eigenvalIter), imag(eigenvalIter)];
% name = sprintf('output/spectrumIter_%s_p%i_k%g_h%g.csv', benchmark, degree, k, h);
% writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');
% 
% disp(['    Min e.v. (Iter)     ' num2str(min(abs(eigenvalIter)))]);
% disp(['    Max e.v. (Iter)     ' num2str(max(abs(eigenvalIter)))]);
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecIter))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecIter))]);
% disp('---------------------------------------------------------');

% figure(1);
% hold off
% scatter(real(eigenvalIter),imag(eigenvalIter),'DisplayName','Eigenvalues');
% hold on
% %plot(fovals(sysA.matS,100),'-b','DisplayName','Numerical range');
% grid on; box on;
% title(['Benchmark "' benchmark '" — k=' num2str(k/pi) 'pi — h=' num2str(degree) ' — h=' num2str(h)]);

end
