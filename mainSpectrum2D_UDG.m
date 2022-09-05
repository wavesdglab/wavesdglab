% close all;
clear all;

headers2D;
global k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE
% benchmark = 'open'; degree = 1; k = 10*pi; h = 1/32;
% benchmark = 'open'; degree = 3; k = 10*pi; h = 1/16;

% BENCH CAVITY
% benchmark = 'cavity'; degree = 1; k = (3+1/8)*sqrt(2)*pi; h = 1/32;
% benchmark = 'cavity'; degree = 3; k = (5+1/8)*sqrt(2)*pi; h = 1/8;

% BENCH WAVEGUIDE
% benchmark = 'waveguide'; degree = 1; k = 2*pi; h = 1/16;
benchmark = 'waveguide'; degree = 3; k = 6*pi; h = 1/8;

tau = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mesh = benchmark2D(benchmark,h);
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

disp(['---------------------------------------------------------']);
disp(['Method UDG (' benchmark ')']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['    tau                 ' num2str(tau)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum2D_UDG(mesh, dofm, tau);
% [errorL2] = computeNormError2D_DG(mesh, dofm, solA);
% 
% [solP, sysP] = computeSolProjL2_2D_DG(mesh, dofm);
% [errorProjL2] = computeNormError2D_DG(mesh, dofm, solP);
% 
% [solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solA);
% [errorPostL2] = computeNormError2D_DG(mesh, dofmPost, solApost);
% 
% [solPpost, sysPpost] = computeSolProjL2_2D_DG(mesh, dofmPost);
% [errorProjPostL2] = computeNormError2D_DG(mesh, dofmPost, solPpost);
% 
% disp(['    L2-Error (numSol)   ' num2str(errorL2, '%1.2e')]);
% disp(['    L2-Error (projSol)  ' num2str(errorProjL2, '%1.2e')]);
% disp(['    L2-Error (numPost)  ' num2str(errorPostL2, '%1.2e')]);
% disp(['    L2-Error (projPost) ' num2str(errorProjPostL2, '%1.2e')]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% writeFieldDG(dofm, mesh, solP, "mySol.pos", "mySol");
% system('gmsh mySol.pos');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('---------------------------------------------------------');
% disp(['Spectrum matA:']);
% disp(['    Size                ' num2str(size(sysA.matA,1))]);
% 
% [eigenvecA,eigenvalA] = eigs(sysA.matA,size(sysA.matA,1));
% eigenvalA = diag(eigenvalA);
% 
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecA))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecA))]);
% disp(['    Cond(A)             ' num2str(condest(sysA.matA))]);
% disp(['    Min eigenval        ' num2str(min(eigenvalA))]);
% 
% figure;
% scatter(real(eigenvalA),imag(eigenvalA),'DisplayName','Eigenvalues');
% %hold on
% %plot(fovals(sysA.matA,100),'-b','DisplayName','Numerical range');
% grid on; box on;
% title(['Eigenvalues matA  —  Min real eigenvalues: ' num2str(min(real(eigenvalA)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 1;
matS = sysA.matS;
matI = sparse(1:size(sysA.matS,1),1:size(sysA.matS,2),1);
matIter = matI - sysA.matS;
mat = (1-alpha)*matI + alpha*matIter;

disp('---------------------------------------------------------');
disp(['Spectrum matS:']);
disp(['    Size                ' num2str(size(mat,1))]);

[eigenvecS,eigenvalS]      = eigs(mat,size(mat,1));
eigenvalS                  = diag(eigenvalS);

% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecS))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecS))]);
% disp(['    Cond(S)             ' num2str(condest(matS), '%1.2e')]);
% disp(['    Min eigenval        ' num2str(min(abs(eigenvalS)))]);
% disp(['    Max eigenval        ' num2str(max(abs(eigenvalS)))]);

figure;
scatter(real(eigenvalS),imag(eigenvalS),'DisplayName','Eigenvalues');
%hold on
%plot(fovals(sysA.matS,100),'-b','DisplayName','Numerical range');
grid on; box on;
title(['Benchmark "' benchmark '" — k=' num2str(k/pi) 'pi — h=' num2str(degree) ' — h=' num2str(h)]);

rezu1 = ["real", "imag"];
rezu2 = [real(eigenvalS), imag(eigenvalS)];
name = sprintf('output/spectrumUDG_%s_P%i_k%g_h%g_tau%g+%gi.csv', benchmark, degree, k, h, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(['---------------------------------------------------------']);
% disp(['Spectrum matPhy:']);
% disp(['    Size                ' num2str(size(sysA.matPhy,1))]);
% 
% [eigenvecP,eigenvalP]      = eigs(sysA.matPhy,size(sysA.matPhy,1));
% eigenvalP                  = diag(eigenvalP);
% 
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecP))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecP))]);
% disp(['    Cond(P)             ' num2str(condest(sysA.matPhy), '%1.2e')]);
% disp(['    Min eigenval        ' num2str(min(eigenvalP))]);
% 
% figure(3);
% hold off;
% scatter(real(eigenvalP),imag(eigenvalP),'DisplayName','Eigenvalues');
% %hold on
% %plot(fovals(sysA.matP,100),'-b','DisplayName','Numerical range');
% grid on; box on;
% title(['Eigenvalues matP  —  Min real eigenvalues: ' num2str(min(real(eigenvalP)))]);