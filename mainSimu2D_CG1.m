%close all;
clear all;

headers2D;

% Define parameters
global k BCWest BCNorth BCEast BCSouth 
k = 30; %2*pi;
h = 0.1; %1/16/2;
degree = 1;

resTol = 1e-4;
BCWest  = 'NEU';
BCNorth = 'NEU';
BCEast  = 'NEU';
BCSouth = 'NEU';

% Build mesh and dofManager
system(['gmsh -2 mesh.geo -v 0 -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh('mesh.msh');
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[solRef, matA, rhsA]     = computeSolNum2D_CG1(mesh, dofm);
solAna                   = computeSolAna2D_CG(mesh);
[errorL2, errorH1]       = computeNormError2D_CG(mesh, dofm, solRef, solAna);
% 
% [solA,~,~,iterA]         = gmres(matA,rhsA,size(matA,1),resTol,size(matA,1));
% errorL2IterA             = computeNormError2D_CG(mesh, dofm, solA, solRef);
% [solA,~,~,iterBiCGStabA] = bicgstab(matA,rhsA,resTol,100*size(matA,1));
% errorL2BiCGStabA         = computeNormError2D_CG(mesh, dofm, solA, solRef);
% [solA,~,~,iterCGNA]      = conjgradn(matA,rhsA,resTol,size(matA,1));
% errorL2CGNA              = computeNormError2D_CG(mesh, dofm, solA, solRef);
% [solA,~,~,iterJacobiA]   = jacobi(matA,rhsA,resTol,10*size(matA,1),0.5);
% errorL2JacobiA           = computeNormError2D_CG(mesh, dofm, solA, solRef);

% [eigenvecA,eigenvalA] = eigs(matA,size(matA,1));
% eigenvalA = diag(eigenvalA);
% [eigenvecAA,eigenvalAA] = eigs(matA'*matA,size(matA,1));
% eigenvalAA = diag(eigenvalAA);

disp(['Method CG-1']);
disp(['---------------------------------------------------------']);
disp(['    L2-Error           ' num2str(errorL2)]);
disp(['    H1-Error           ' num2str(errorH1)]);
disp(['---------------------------------------------------------']);
disp(['A : Size               ' num2str(size(matA,1))]);
% disp(['    Rank(eigenvectors) ' num2str(rank(eigenvecA))]);
% disp(['    Cond(eigenvectors) ' num2str(cond(eigenvecA))]);
% disp(['    Cond(A)            ' num2str(condest(matA))]);
% disp(['    - - - - - - - - - - - - - - - - - - - - -']);
% disp(['    IterGmres          ' num2str(iterA(2))]);
% disp(['    Final L2-Error     ' num2str(errorL2IterA)]);
% disp(['    IterBiCGS          ' num2str(iterBiCGStabA)]);
% disp(['    Final L2-Error     ' num2str(errorL2BiCGStabA)]);
% disp(['    IterCGN            ' num2str(iterCGNA)]);
% disp(['    Final L2-Error     ' num2str(errorL2CGNA)]);
% disp(['    IterRelax          ' num2str(iterJacobiA)]);
% disp(['    Final L2-Error     ' num2str(errorL2JacobiA)]);
disp(['---------------------------------------------------------']);

% disp(['\text{CG-1} & & ' ...
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

% figure(1);
% subplot(1,3,1);
% postProVizuCG(mesh, real(solRef), 'Exact solution');
% subplot(1,3,2);
% postProVizuCG(mesh, real(solA), 'Numerical solution');
% subplot(1,3,3);
% postProVizuCG(mesh, real(solA-solRef), 'Error');