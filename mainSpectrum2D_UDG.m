%close all;
clear all;

headers2D;
global k;

tau = 1;
prec = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE
benchmark = 'open'; degree = 3; k = 15*pi; h = 1/16;
run(benchmark,degree,h,tau,prec);
benchmark = 'open'; degree = 3; k = 15*pi; h = 1/34;
run(benchmark,degree,h,tau,prec);
benchmark = 'open'; degree = 3; k = 30*pi; h = 1/34;
run(benchmark,degree,h,tau,prec);

% BENCH CAVITY
benchmark = 'cavity'; degree = 3; k = (7+1/10)*sqrt(2)*pi; h = 1/10;
run(benchmark,degree,h,tau,prec);
benchmark = 'cavity'; degree = 3; k = (7+1/10)*sqrt(2)*pi; h = 1/15;
run(benchmark,degree,h,tau,prec);
benchmark = 'cavity'; degree = 3; k = (7+1/100)*sqrt(2)*pi; h = 1/15;
run(benchmark,degree,h,tau,prec);

% BENCH WAVEGUIDE
benchmark = 'waveguide'; degree = 3; k = 6*pi; h = 1/8;
run(benchmark,degree,h,tau,prec);
benchmark = 'waveguide'; degree = 3; k = 6*pi; h = 1/17;
run(benchmark,degree,h,tau,prec);
benchmark = 'waveguide'; degree = 3; k = 12*pi; h = 1/17;
run(benchmark,degree,h,tau,prec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,h,tau,prec)
global k;

mesh = benchmark2D(benchmark,h);
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

[solA, sysA, condLoc] = computeSolNum2D_UDG(mesh, dofm, tau, prec);
[normErr, normErrU, normErrV, normSol, normSolU, normSolV] = computeNormError2D_DG(mesh, dofm, solA);

[solP] = computeSolProjL2_2D_DG(mesh, dofm);
[normProjErr, normProjErrU, normProjErrV] = computeNormError2D_DG(mesh, dofm, solP);

disp(['    L2-Norm Sol       ' num2str(normSol, '%1.2e') '  ' num2str(normSolU, '%1.2e') '  ' num2str(normSolV, '%1.2e')]);
disp(['    L2-Norm ErrorSol  ' num2str(normErr, '%1.2e') '  ' num2str(normErrU, '%1.2e') '  ' num2str(normErrV, '%1.2e')]);
disp(['    L2-Norm ErrorProj ' num2str(normProjErr, '%1.2e') '  ' num2str(normProjErrU, '%1.2e') '  ' num2str(normProjErrV, '%1.2e')]);
disp('---------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% writeFieldDG(dofm, mesh, solP, "mySol.pos", "mySol");
% system('gmsh mySol.pos');

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
name = sprintf('output/statsCHDG_%s.csv', benchmark);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alpha = 1;
% matI = sparse(1:size(sysA.matS,1),1:size(sysA.matS,2),1);
% matIter = (1-alpha)*matI + alpha*(matI-sysA.matS);
% 
% [eigenvecGG,eigenvalGG] = eigs(sysA.matGG,size(sysA.matGG,1));
% 
% precL = inv(sysA.matGG);
% precR = sparse(1:size(sysA.matGG,1),1:size(sysA.matGG,1),1);
% 
% precL = eigenvalGG\eigenvecGG';
% precR = eigenvecGG;
% 
% matIter = precL*sysA.matS*precR;
% 
% eigenvalI = eigs(matIter,size(matIter,1));
% 
% [eigenvecI,eigenvalI] = eigs(matIter,size(matIter,1));
% eigenvalI            = diag(eigenvalI);
% 
% disp(['    Spectral radius     ' num2str(max(abs(eigenvalI)))]);


% disp('---------------------------------------------------------');
% disp(['Reduced matrix:']);
% 
% [eigenvecS,eigenvalS] = eigs(sysA.matS,size(sysA.matS,1));
% eigenvalS             = diag(eigenvalS);
% 
% disp(['    Min e.v. (Iter)     ' num2str(min(abs(eigenvalS)))]);
% disp(['    Max e.v. (Iter)     ' num2str(max(abs(eigenvalS)))]);
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecS))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecS))]);
% 
% disp('---------------------------------------------------------');
% disp(['Conditionning:']);
% 
% condPhy = condest(sysA.matPhy);
% disp(['    Condest(Phy)        ' num2str(condPhy, '%1.2e')]);
% condRed = condest(sysA.matS);
% disp(['    Condest(Red)        ' num2str(condRed, '%1.2e')]);
% condLocMin = min(condLoc);
% disp(['    CondMin(Loc)        ' num2str(condLocMin, '%1.2e')]);
% condLocMax = max(condLoc);
% disp(['    CondMax(Loc)        ' num2str(condLocMax, '%1.2e')]);
% 
% disp('---------------------------------------------------------');

% figure(3);
% scatter(real(eigenvalI),imag(eigenvalI),'DisplayName','Eigenvalues');
% hold on
% plot(fovals(sysA.matS,100),'-b','DisplayName','Numerical range');
% grid on; box on;
% title(['Benchmark "' benchmark '" — k=' num2str(k/pi) 'pi — h=' num2str(degree) ' — h=' num2str(h)]);

% rezu1 = ["real", "imag"];
% rezu2 = [real(eigenvalI), imag(eigenvalI)];
% name = sprintf('output/spectrumI_%s_P%i_k%g_h%g_tau%g+%gi.csv', benchmark, degree, k, h, real(tau), imag(tau));
% writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');
% 
% rezu1 = ["real", "imag"];
% rezu2 = [real(eigenvalS), imag(eigenvalS)];
% name = sprintf('output/spectrumUDG_%s_P%i_k%g_h%g_tau%g+%gi.csv', benchmark, degree, k, h, real(tau), imag(tau));
% writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

end
