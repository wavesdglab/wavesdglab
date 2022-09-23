% close all;
clear all;

headers2D;
global k

tau = 1;
theta = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE
benchmark = 'open'; degree = 3; k = 15*pi; h = 1/16;
run(benchmark,degree,h,tau,theta)
% benchmark = 'open'; degree = 3; k = 15*pi; h = 1/16/2;
% run(benchmark,degree,h,tau,theta)
% benchmark = 'open'; degree = 3; k = 15*pi*2; h = 1/16/2;
% run(benchmark,degree,h,tau,theta)

% BENCH CAVITY
benchmark = 'cavity'; degree = 3; k = (5+1/8)*sqrt(2)*pi; h = 1/8;
run(benchmark,degree,h,tau,theta)
% benchmark = 'cavity'; degree = 3; k = (5+1/8)*sqrt(2)*pi; h = 1/8/2;
% run(benchmark,degree,h,tau,theta)
% benchmark = 'cavity'; degree = 3; k = (10+1/8)*sqrt(2)*pi; h = 1/8/2;
% run(benchmark,degree,h,tau,theta)

% BENCH WAVEGUIDE
benchmark = 'waveguide'; degree = 3; k = 6*pi; h = 1/8;
run(benchmark,degree,h,tau,theta)
% benchmark = 'waveguide'; degree = 3; k = 6*pi; h = 1/8/2;
% run(benchmark,degree,h,tau,theta)
% benchmark = 'waveguide'; degree = 3; k = 6*pi*2; h = 1/8/2;
% run(benchmark,degree,h,tau,theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,h,tau,theta)
global k;

mesh = benchmark2D(benchmark,h);
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

disp(['=========================================================']);
disp(['Method DG (' benchmark ')']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['    tau                 ' num2str(tau)]);
disp(['    theta               ' num2str(theta)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_DG(mesh, dofm, tau, theta);
[errorL2] = computeNormError2D_DG(mesh, dofm, solA);

[solP, ~] = computeSolProjL2_2D_DG(mesh, dofm);
[errorProjL2] = computeNormError2D_DG(mesh, dofm, solP);

[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solA);
[errorPostL2] = computeNormError2D_DG(mesh, dofmPost, solApost);

[solPpost, ~] = computeSolProjL2_2D_DG(mesh, dofmPost);
[errorProjPostL2] = computeNormError2D_DG(mesh, dofmPost, solPpost);

disp(['    L2-Error (numSol)   ' num2str(errorL2, '%1.2e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2, '%1.2e')]);
disp(['    L2-Error (numPost)  ' num2str(errorPostL2, '%1.2e')]);
disp(['    L2-Error (projPost) ' num2str(errorProjPostL2, '%1.2e')]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% writeFieldDG(dofm, mesh, solP, "mySol.pos", "mySol");
% system('gmsh mySol.pos');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------------------------------------------');
disp(['    Size                ' num2str(size(sysA.matA,1))]);
disp(['    nnz(S)              ' num2str(nnz(sysA.matA))]);

%[eigenvecA,eigenvalA] = eigs(sysA.matA,size(sysA.matA,1));
%eigenvalA             = diag(eigenvalA);

%disp(['    Min e.v. (Iter)     ' num2str(min(abs(eigenvalA)))]);
%disp(['    Max e.v. (Iter)     ' num2str(max(abs(eigenvalA)))]);
%disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecA))]);
%disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecA))]);
disp(['    Condest(A)          ' num2str(condest(sysA.matA), '%1.2e')]);

disp('---------------------------------------------------------');

% figure(2);
% hold off
% scatter(real(eigenvalS),imag(eigenvalS),'DisplayName','Eigenvalues');
% %plot(fovals(sysA.matS,100),'-b','DisplayName','Numerical range');
% grid on; box on;
% title(['Benchmark "' benchmark '" — k=' num2str(k/pi) 'pi — h=' num2str(degree) ' — h=' num2str(h)]);

% rezu1 = ["real", "imag"];
% rezu2 = [real(eigenvalA), imag(eigenvalA)];
% name = sprintf('output/spectrumHDG_%s_P%i_k%g_h%g_tau%g+%gi.csv', benchmark, degree, k, h, real(tau,theta), imag(tau,theta));
% writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

end
