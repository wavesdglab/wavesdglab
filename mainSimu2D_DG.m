%close all;
clear all;

headers2D;

% Define parameters
global k BCWest BCNorth BCEast BCSouth 
k = 3;
h = 0.1;
degree = 1;
tau = 1;
resTol = 1e-4;
BCWest  = 'ABC';
BCNorth = 'ABC';
BCEast  = 'ABC';
BCSouth = 'ABC';

% Build mesh and dofManager
system(['gmsh -2 mesh.geo -v 0 -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh('mesh.msh');
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[solRef, matA, rhsA]     = computeSolNum2D_DG1(mesh, dofm, tau);
solAna                   = computeSolAna2D_DG(mesh);
[errorL2, errorH1]       = computeNormError2D_DG(mesh, dofm, solRef, solAna);

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

% [eigenvecA,eigenvalA] = eigs(matA,size(matA,1));
% eigenvalA = diag(eigenvalA);
% [eigenvecAA,eigenvalAA] = eigs(matA'*matA,size(matA,1));
% eigenvalAA = diag(eigenvalAA);

disp(['Method DG-1']);
disp(['---------------------------------------------------------']);
disp(['    L2-Error           ' num2str(errorL2)]);
disp(['    H1-Error           ' num2str(errorH1)]);
disp(['---------------------------------------------------------']);
% disp(['A : Size               ' num2str(size(matA,1))]);
% disp(['    Rank(eigenvectors) ' num2str(rank(eigenvecA))]);
% disp(['    Cond(eigenvectors) ' num2str(cond(eigenvecA))]);
% disp(['    Cond(A)            ' num2str(condest(matA))]);
% disp(['    Cond(AA)           ' num2str(condest(matA'*matA))]);
% disp(['    Min eigenval AA    ' num2str(min(eigenvalAA))]);
% disp(['    - - - - - - - - - - - - - - - - - - - - -']);
% disp(['    IterGmres          ' num2str(iterA(2))]);
% disp(['    Final L2-Error     ' num2str(errorL2IterA)]);
% disp(['    IterBiCGS          ' num2str(iterBiCGStabA)]);
% disp(['    Final L2-Error     ' num2str(errorL2BiCGStabA)]);
% disp(['    IterCGN            ' num2str(iterCGNA)]);
% disp(['    Final L2-Error     ' num2str(errorL2CGNA)]);
% disp(['    IterRelax          ' num2str(iterJacobiA)]);
% disp(['    Final L2-Error     ' num2str(errorL2JacobiA)]);
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

figure(1);
subplot(1,3,1);
hold off
postProVizu2D_DG(mesh,real(solAna), 'Analytic solution');
subplot(1,3,2);
hold off
postProVizu2D_DG(mesh,real(solRef), 'Numerical solution');
subplot(1,3,3);
hold off
postProVizu2D_DG(mesh,real(solRef-solAna), 'Error');
