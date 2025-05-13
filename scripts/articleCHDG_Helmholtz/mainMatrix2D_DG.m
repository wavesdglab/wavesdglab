%close all;
clear;

global k h

degree = 3;
tau = 1;
theta = 1;
PREC = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE
benchmark = 'open'; k = 15*pi; h = 1/16;
run(benchmark,degree,tau,theta,PREC);
benchmark = 'open'; degree = 3; k = 15*pi; h = 1/34;
run(benchmark,degree,tau,theta,PREC);
benchmark = 'open'; degree = 3; k = 30*pi; h = 1/34;
run(benchmark,degree,tau,theta,PREC);

% BENCH CAVITY
benchmark = 'cavity'; k = (7+1/10)*sqrt(2)*pi; h = 1/10;
run(benchmark,degree,tau,theta,PREC);
benchmark = 'cavity'; degree = 3; k = (7+1/10)*sqrt(2)*pi; h = 1/15;
run(benchmark,degree,tau,theta,PREC);
benchmark = 'cavity'; degree = 3; k = (7+1/100)*sqrt(2)*pi; h = 1/15;
run(benchmark,degree,tau,theta,PREC);

% BENCH WAVEGUIDE
benchmark = 'waveguide'; k = 6*pi; h = 1/8;
run(benchmark,degree,tau,theta,PREC);
benchmark = 'waveguide'; degree = 3; k = 6*pi; h = 1/17;
run(benchmark,degree,tau,theta,PREC);
benchmark = 'waveguide'; degree = 3; k = 12*pi; h = 1/17;
run(benchmark,degree,tau,theta,PREC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,tau,theta,PREC)
global k h

mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method DG (' benchmark ')']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_DG(mesh, dofm, tau, theta, PREC);
[normErr, normSol] = computeNormError2D_DG(mesh, dofm, solA);

solP = computeSolProjL2_2D_DG(mesh, dofm);
normProjErr = computeNormError2D_DG(mesh, dofm, solP);

disp(['    L2-Norm Sol       ' num2str(normSol, '%1.2e')]);
disp(['    L2-Norm ErrorSol  ' num2str(normErr, '%1.2e')]);
disp(['    L2-Norm ErrorProj ' num2str(normProjErr, '%1.2e')]);
disp(['---------------------------------------------------------']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sizeA = size(sysA.matA,1);
disp(['    Size                ' num2str(sizeA)]);
nnzA = nnz(sysA.matA);
disp(['    nnz(S)              ' num2str(nnzA)]);
condA = condest(sysA.matA);
disp(['    Condest(A)          ' num2str(condA, '%1.2e')]);
disp(['---------------------------------------------------------']);

rezu1 = ["degree" "k" "h" "real(tau)" "imag(tau)" "normErr" "normProjErr" "normSol" "sizeA" "nnzA" "condA"];
rezu2 = [degree, k, h, real(tau), imag(tau), normErr, normProjErr, normSol, sizeA, nnzA, condA];
name = sprintf('output/statsDG_%s_p%i_k%g_h%g.csv', benchmark, degree, k, h);
writematrix([rezu1 ; rezu2],name,'Delimiter','semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [eigenvecA,eigenvalA] = eigs(sysA.matA,size(sysA.matA,1));
% eigenvalA             = diag(eigenvalA);
% 
% disp(['    Min e.v. (Iter)     ' num2str(min(abs(eigenvalA)))]);
% disp(['    Max e.v. (Iter)     ' num2str(max(abs(eigenvalA)))]);
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecA))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecA))]);
% 
% figure(2);
% hold off
% scatter(real(eigenvalS),imag(eigenvalS),'DisplayName','Eigenvalues');
% %plot(fovals(sysA.matS,100),'-b','DisplayName','Numerical range');
% grid on; box on;
% title(['Benchmark "' benchmark '" — k=' num2str(k/pi) 'pi — h=' num2str(degree) ' — h=' num2str(h)]);
% 
% rezu1 = ["real", "imag"];
% rezu2 = [real(eigenvalA), imag(eigenvalA)];
% name = sprintf('output/spectrumHDG_%s_P%i_k%g_h%g_tau%g+%gi.csv', benchmark, degree, k, h, real(tau,theta), imag(tau,theta));
% writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

end
