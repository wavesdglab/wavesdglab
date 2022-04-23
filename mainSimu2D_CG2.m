%close all;
clear all;

headers2D;

% Define parameters
global k
k = 50;       % 2*pi;
h = 0.1;      % 1/16/2;
degree = 4;

resTol = 1e-4;
benchmark2D('open');

% Build mesh and dofManager
system(['gmsh -2 mesh.geo -v 0 -clmax ' num2str(h) ' -clmin ' num2str(h)]);
mesh = readMesh('mesh.msh');
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['---------------------------------------------------------']);
disp(['Method CG-2']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

[solA, matA, rhsA] = computeSolNum2D_CG2(mesh, dofm);
[errorL2, errorH1] = computeNormError2D_CG(mesh, dofm, solA);

[solP, matP, rhsP] = computeSolProjL2_2D_CG(mesh, dofm);
[errorProjL2, errorProjH1] = computeNormError2D_CG(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(errorL2)]);
%disp(['    H1-Error (numSol)   ' num2str(errorH1)]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2)]);
%disp(['    H1-Error (projSol)  ' num2str(errorProjH1)]);
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

% numDofTRIred = mesh.numVer * dofm.numDofPerVer + mesh.numEdg * dofm.numDofPerEdg;
% dofG = 1:numDofTRIred;
% dofI = (numDofTRIred+1):dofm.numDofTRI;
% matGG = matA(dofG,dofG);
% matGI = matA(dofG,dofI);
% matIG = matA(dofI,dofG);
% matII = matA(dofI,dofI);
% rhsG = rhsA(dofG);
% rhsI = rhsA(dofI);
% matS = matGG - matGI*(matII\matIG);
% rhsS = rhsG - matGI*(matII\rhsI);
% solG = matS\rhsS;
% solI = matII\(rhsI-matIG*solG);
% solRedu = [ solG ; solI ];
% 
% [eigenvecS,eigenvalS]      = eigs(matS,size(matS,1));
% eigenvalS                  = diag(eigenvalS);
% 
% disp(['S : Size                ' num2str(size(matS,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rank(eigenvecS))]);
% disp(['    Cond(eigenvectors)  ' num2str(cond(eigenvecS))]);
% disp(['    Cond(S)             ' num2str(condest(matS))]);
% disp(['---------------------------------------------------------']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(['--- CALL gmres']);
% [solAiter,~,~,iterA]         = gmres(matA,rhsA,size(matA,1),resTol,size(matA,1));
% errorL2IterA                 = computeNormError2D_CG(mesh, dofm, solAiter, solA);
% disp(['--- CALL bicgstab']);
% [solAiter,~,~,iterBiCGStabA] = bicgstab(matA,rhsA,resTol,100*size(matA,1));
% errorL2BiCGStabA             = computeNormError2D_CG(mesh, dofm, solAiter, solA);
% disp(['--- CALL conjgradn']);
% [solAiter,~,~,iterCGNA]      = conjgradn(matA,rhsA,resTol,size(matA,1));
% errorL2CGNA                  = computeNormError2D_CG(mesh, dofm, solAiter, solA);
% disp(['--- CALL jacobi']);
% [solAiter,~,~,iterJacobiA]   = jacobi(matA,rhsA,resTol,size(matA,1),0.5);
% errorL2JacobiA               = computeNormError2D_CG(mesh, dofm, solAiter, solA);
% disp(['--- CALL eigenvalues']);
% 
% disp(['---------------------------------------------------------']);
% disp(['    IterGmres           ' num2str(iterA(2))]);
% disp(['    Final L2-Error      ' num2str(errorL2IterA)]);
% disp(['    IterBiCGS           ' num2str(iterBiCGStabA)]);
% disp(['    Final L2-Error      ' num2str(errorL2BiCGStabA)]);
% disp(['    IterCGN             ' num2str(iterCGNA)]);
% disp(['    Final L2-Error      ' num2str(errorL2CGNA)]);
% disp(['    IterRelax           ' num2str(iterJacobiA)]);
% disp(['    Final L2-Error      ' num2str(errorL2JacobiA)]);
% disp(['---------------------------------------------------------']);

% disp(['\text{CG-2} & & ' ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% writeFieldCG(dofm, mesh, solA, "mySol.pos", "mySol");
% system('gmsh mySol.pos');

% figure(1);
% subplot(1,3,1);
% hold off;
% postProVizu2D_CG(mesh, real(solP), 'Projected solution');
% subplot(1,3,2);
% hold off;
% postProVizu2D_CG(mesh, real(solA), 'Numerical solution');
% subplot(1,3,3);
% hold off;
% postProVizu2D_CG(mesh, real(solA-solP), 'Error');
