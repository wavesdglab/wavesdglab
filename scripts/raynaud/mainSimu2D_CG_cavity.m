clear all;
%close all;

global k h

% Setup benchmark and parameters
benchmark = 'cavity';
% k = 1.1*sqrt(2)*pi; h = 1/10;
k = 1.00001*sqrt(10)*pi; h = 1/16;
tol = 1e-10; maxit = 2000; itout = 5;
degree = 3;
PREC = 0;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CG - Benchmark "' benchmark '"']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

% Compute numerical solution/error
[solA, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);
% errorL2 = computeNormError2D_CG(mesh, dofm, solA);
% disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);

% Compute projection solution/error
% solP = computeSolProjL2_2D_CG(mesh, dofm);
% errorProjL2 = computeNormError2D_CG(mesh, dofm, solP);
% disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
% disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

% writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
% writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
% global PML_HIDE; PML_HIDE = 1;
% writeField2D(dofm, mesh, solA-solP, 'output/errNum.pos', "errNum");


% -------------------------------------------------------------------------
% Compute eigenvalues/eigenvectors
% -------------------------------------------------------------------------

mat = sysA.matA;
[eigenvec, eigenval] = eigs(mat,20,'smallestabs');
eigenval = diag(eigenval);
% rankEigenVec = rank(eigenvec);
% condEigenVec = cond(eigenvec);

csvwrite('output/eigenvec.csv',eigenvec);

% 
% disp(['    Size                ' num2str(size(mat,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rankEigenVec)]);
% disp(['    Cond(eigenvectors)  ' num2str(condEigenVec)]);
% 
% % Plot spectrum
% figure;
% hold off; scatter(real(eigenval),imag(eigenval));
% % hold on; plot(fovals(mat,100));
% grid on; box on;
% 
% % Compute condition number
% condestMat = condest(mat);
% 
% disp(['    Cond(mat)           ' num2str(condestMat,'%1.2e')]);
% disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Compute iterative solution
% -------------------------------------------------------------------------

solver = 'GMRES_dev';
[resVec, ~, i_cv, ~, ~, X, hrv, ~, dist] = solverGMRES_dev(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);

resVec = resVec(1:i_cv/itout+1);
dist = dist(2:i_cv/itout+1);
dist = dist';
hrv = hrv(2:i_cv/itout+1,:);

csvwrite('output/resVec.csv',resVec);
csvwrite('output/dist.csv',dist);
csvwrite('output/hrv.csv',hrv);



% x=0:itout:i_cv;
% solGMRES = X(:,end);
% writeField2D(dofm, mesh, solGMRES, 'output/solGMRES.pos', "solGMRES");
% system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos output/solGMRES.pos&');
% 
% figure;
% 
% semilogy(x,resVec,'k-o');
% hold off
% legend('Relative residual');
% box on;
% grid on;
% legend('Location','southwest');
% title(['CG - ' benchmark ' - ' solver ' - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)])
% xlim([0 i_cv]);
% ylim([tol 1]);
% xlabel('Iteration');
% ylabel('Value');