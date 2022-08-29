%close all;
clear all;

headers2D;
global k

% hPower = [-1 -1.5 -2 -2.5 -3 -3.5];
% hList = 2.^hPower;
% for ITER = 1:size(hPower,2)
% h = hList(ITER);

% FREE SPACE
k = 10*pi;
h = 1/8;
degree = 3;
mesh = benchmark2D('open',h);

% WAVEGUIDE
% k = 6*pi;
% h = 1/8;
% degree = 3;
% mesh = benchmark2D('waveguide',h);

% CAVITY
% k = 5.125*sqrt(2)*pi;
% h = 1/8;
% degree = 3;
% mesh = benchmark2D('cavity',h);

% Define parameters
resTol = 1e-4;
tau = 1;
theta = 0;

% Build mesh and dofManager
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method DG-1']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

[solA, matA, rhsA] = computeSolNum2D_DG1(mesh, dofm, tau, theta);
[errorL2, errorH1] = computeNormError2D_DG(mesh, dofm, solA);

[solP, sysP] = computeSolProjL2_2D(mesh, dofm);
[errorProjL2, errorProjH1] = computeNormError2D_DG(mesh, dofm, solP);

[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solA);
[errorPostL2, errorPostH1] = computeNormError2D_DG(mesh, dofmPost, solApost);

[solPpost, sysPpost] = computeSolProjL2_2D(mesh, dofmPost);
[errorProjPostL2, errorProjPostH1] = computeNormError2D_DG(mesh, dofmPost, solPpost);

% disp([num2str(k) ' ' num2str(h) ' ' num2str(degree) ' ' num2str(Dlambda) ' ' num2str(errorL2) ' ' num2str(errorProjL2) ' ' num2str(errorPostL2) ' ' num2str(errorProjPostL2)]);

disp(['    L2-Error (numSol)   ' num2str(errorL2)]);
% disp(['    H1-Error (numSol)   ' num2str(errorH1)]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2)]);
% disp(['    H1-Error (projSol)  ' num2str(errorProjH1)]);
disp(['---------------------------------------------------------']);

% [eigenvecA,eigenvalA] = eigs(matA,size(matA,1));
% eigenvalA = diag(eigenvalA);
% 
% disp(['A : Size                ' num2str(size(matA,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecA))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecA))]);
% disp(['    Cond(A)             ' num2str(condest(matA))]);
% disp(['---------------------------------------------------------']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

writeFieldDG(dofm, mesh, solA, "mySol.pos", "mySol");
system('gmsh mySol.pos');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('Solver  : gmres A\n');
% [solA,~,~,iterA]         = gmres(matA,rhsA,size(matA,1),resTol,size(matA,1));
% errorL2IterA             = computeNormError2D_DG(mesh, dofm, solA, solRef);
% fprintf('Solver  : bicgstab A\n');
% [solA,~,~,iterBiCGStabA] = bicgstab(matA,rhsA,resTol,size(matA,1));
% errorL2BiCGStabA         = computeNormError2D_DG(mesh, dofm, solA, solRef);
% fprintf('Solver  : conjgradn A\n');
% [solA,~,~,iterCGNA]      = conjgradn(matA,rhsA,resTol,size(matA,1));
% errorL2CGNA              = computeNormError2D_DG(mesh, dofm, solA, solRef);
% fprintf('Solver  : relaxation A\n');
% [solA,~,~,iterJacobiA]   = jacobi(matA,rhsA,resTol,size(matA,1),0.5);
% errorL2JacobiA           = computeNormError2D_DG(mesh, dofm, solA, solRef);
% %iterJacobiA = NaN;
% %errorL2JacobiA = NaN;

% fprintf('Solver  : eigenval A\n');
% [eigenvecA,eigenvalA] = eigs(matA,size(matA,1));
% eigenvalA = diag(eigenvalA);
% [eigenvecAA,eigenvalAA] = eigs(matA'*matA,size(matA,1));
% eigenvalAA = diag(eigenvalAA);

% figure;
% hold off
% scatter(real(eigenvalA),imag(eigenvalA),'b','DisplayName','Eigenvalues');
% hold on
% %plot(fovals(matS,100),'-b','DisplayName','Numerical range');
% grid on; box on;
% title(['Eigenvalues : ' BCWest ' + ' BCNorth ' + ' BCEast ' + ' BCSouth]);
% legend();
% %axis([-0.1 1.1 -1.5 1.2]);
% %axis([-0.05 0.5 -1.5 1.2]);

% disp(['A : Size                ' num2str(size(matA,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecA))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecA))]);
% disp(['    Cond(A)             ' num2str(condest(matA))]);
% disp(['    Cond(AA)            ' num2str(condest(matA'*matA))]);
% disp(['    Min eigenval AA     ' num2str(min(eigenvalAA))]);
% disp(['    - - - - - - - - - - - - - - - - - - - - -']);
% disp(['    IterGmres           ' num2str(iterA(2))]);
% disp(['    Final L2-Error      ' num2str(errorL2IterA)]);
% disp(['    IterBiCGS           ' num2str(iterBiCGStabA)]);
% disp(['    Final L2-Error      ' num2str(errorL2BiCGStabA)]);
% disp(['    IterCGN              ' num2str(iterCGNA)]);
% disp(['    Final L2-Error      ' num2str(errorL2CGNA)]);
% disp(['    IterRelax           ' num2str(iterJacobiA)]);
% disp(['    Final L2-Error      ' num2str(errorL2JacobiA)]);
% disp(['---------------------------------------------------------']);

% disp(['\text{DG} & & ' ...
%     num2str(errorL2,'%.1e') ' & ' ...
%     num2str(errorH1,'%.1e') ' & ' ...
%     num2str(size(matA,1)) ' & ' ...
%     num2str(rank(eigenvecA)) ' & ' ...
%     num2str(cond(eigenvecA),'%.1e') ' & ' ...
%     num2str(condest(matA),'%.1e') ' & ' ...
%     num2str(iterA(2)) ' & ' ...
%     num2str(errorL2IterA,'%.1e') ' & ' ...
%     num2str(iterBiCGStabA) ' & ' ...
%     num2str(errorL2BiCGStabA,'%.1e') ' & ' ...
%     num2str(iterCGNA) ' & ' ...
%     num2str(errorL2CGNA,'%.1e') ' & ' ...
%     num2str(iterJacobiA) ' & ' ...
%     num2str(errorL2JacobiA,'%.1e') ' \\'
%     ]);

% end
