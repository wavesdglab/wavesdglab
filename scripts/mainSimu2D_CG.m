clear all;
%close all;

global k

% Setup benchmark and parameters
benchmark = 'cavity2';
switch benchmark
    case 'open'
        k = 15*pi;
        h = 1/16;
        tol = 1e-10; maxit = 1000; itout = 50;
    case 'cavity'
        k = 2*sqrt(2)*pi;
        h = 1/64;
        tol = 1e-10; maxit = 2000; itout = 100;
    case 'cavity2'
        k = sqrt(5+0.5)*pi;
        h = 1/64;
        tol = 1e-10; maxit = 2000; itout = 100;
    case 'waveguide'
        k = 6*pi;
        h = 1/8;
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 1; % P1
PREC = 0; % for preconditioner

% Build mesh and DOF manager
mesh = benchmark2D(benchmark,h);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree); % espace fonctionnel discret

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1); % nb de points par longueur d'onde

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

[solA, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);
errorL2 = computeNormError2D_CG(mesh, dofm, solA);
% 
solP = computeSolProjL2_2D_CG(mesh, dofm);
errorProjL2 = computeNormError2D_CG(mesh, dofm, solP);
% 
% disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
% disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
% disp(['---------------------------------------------------------']);

% disp([num2str(k) ' ' num2str(h) ' ' num2str(degree) ' ' num2str(Dlambda) ' ' num2str(errorL2) ' ' num2str(errorProjL2)]);

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
writeField2D(dofm, mesh, solA-solP, 'output/errNum.pos', "errNum");
system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

% -------------------------------------------------------------------------
% Compute eigenvalues/eigenvectors
% -------------------------------------------------------------------------

% mat = sysA.matA;
% [eigenvec, eigenval] = eigs(mat,size(mat,1));
% eigenval = diag(eigenval);
% rankEigenVec = rank(eigenvec);
% condEigenVec = cond(eigenvec);

% writeField2D(dofm, mesh, eigenvec(:,end), 'output/eigenvec.pos', "eigenvec");
% system('gmsh output/eigenvec.pos');


% disp(['    Size                ' num2str(size(mat,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rankEigenVec)]);
% disp(['    Cond(eigenvectors)  ' num2str(condEigenVec)]);

% Plot spectrum
% figure;
% hold off; scatter(real(eigenval),imag(eigenval));
% % hold on; plot(fovals(mat,100));
% grid on; box on;

% Compute condition number
% condestMat = condest(mat);
% 
% disp(['    Cond(mat)           ' num2str(condestMat,'%1.2e')]);
% disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Compute iterative solution
% -------------------------------------------------------------------------

% solver = 'GMRES';
% switch solver
%     case 'CGNR'
%         [resRedVec, resPhyVec, errorVec] = solverCGNRredu_CG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
%     case 'GMRES'
%         [resVec, errorVec] = solverGMRES(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
%         %[resRedVec, resPhyVec, errorVec] = solverGMRESredu_CG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
% end

% figure;
% hold off
%semilogy(0:itout:maxit,resPhyVec,'r-o','DisplayName','Relative residual (Phy)');
%semilogy(0:itout:maxit,resRedVec,'r-o','DisplayName','Relative residual (Red)');
% hold on
% semilogy(0:itout:maxit,resVec,'b-','DisplayName','Relative residual');
% semilogy(0:itout:maxit,errorVec,'k','DisplayName','Relative L2-error (iterative)');
%plot([0 maxit],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (direct)');
%plot([0 maxit],[errorProjL2 errorProjL2],'k:','DisplayName','Relative L2-error (projection)');
% box on;
% grid on;
% legend('Location','southwest');
% title(['CG - ' benchmark ' - ' solver ' - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)])
% xlim([0 maxit]);
% ylim([0.005 1]);
% xlabel('Iteration');
% ylabel('Value');