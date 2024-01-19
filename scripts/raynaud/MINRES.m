clear all;
%close all;

global k L L_PML R_disk l;
L_PML = 0;
R_disk = 0;
l = 0;


N=10;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);

computeSolNum2D = @computeSolNum2D_CG;

% Setup benchmark and parameters
benchmark = 'cavity';
switch benchmark
    case 'open'
        k = 15*pi;
        h = 1/16;
        tol = 1e-10; maxit = 1000; itout = 50;
    case 'cavity'
        k = 2.001*sqrt(2)*pi;
%         k = 5.877;
        h = 1/8;
        tol = 1e-10; maxit = 2000; itout =1;
        L = 1;
    case 'scatteringPML'
        k = 25;
        h = 0.05;
        tol = 1e-10; maxit = 2000; itout = 50;
        L = 1.1;
        R_disk = 1;
        L_PML = 0.2;
        computeSolNum2D = @computeSolNum2D_PML_CG;
    case 'scattering_rect'
        k = 7.5*pi;
        h = 0.4;
        tol = 1e-7; maxit = 5000; itout = 10;
        L = 0.95;
        L_PML = 0.2;
        l = 0.5;
        computeSolNum2D = @computeSolNum2D_PML_CG;
    case 'waveguide'
        k = 6*pi;
        h = 1/8;
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 3; % P1
PREC = 1; % for preconditioner

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

[~, sysA] = computeSolNum2D(mesh, dofm, PREC);
% errorL2 = computeNormError2D_CG(mesh, dofm, solA);
%
% solP = computeSolProjL2_2D_CG(mesh, dofm);
% errorProjL2 = computeNormError2D_CG(mesh, dofm, solP);
% 
% disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
% disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
% disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Write and vizu solution
% -------------------------------------------------------------------------

% writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
% writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
% writeField2D(dofm, mesh, solA-solP, 'output/errNum.pos', "errNum");
% system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

% -------------------------------------------------------------------------
% Compute eigenvalues/eigenvectors
% -------------------------------------------------------------------------

% mat = sysA.matPinv * sysA.matA;
A = sysA.matP\sysA.matA;
mat = blkdiag(A,-A');
[eigenvec, eigenval] = eigs(mat,size(mat,1));
eigenval = diag(eigenval);
% return;
% rankEigenVec = rank(eigenvec);
% condEigenVec = cond(eigenvec);
% get the largest eigenvalue in absolute value
lambdaMax = max(abs(eigenval));
% get the smallest eigenvalue in absolute value
lambdaMin = min(abs(eigenval));


% disp(['    Size                ' num2str(size(mat,1))]);
% disp(['    Rank(eigenvectors)  ' num2str(rankEigenVec)]);
% disp(['    Cond(eigenvectors)  ' num2str(condEigenVec)]);

% csvwrite(['output/eigA_' num2str(k) '.csv'], eigenval');

% Plot spectrum
% figure;
% hold off; scatter(real(eigenval),imag(eigenval));
% % hold on; plot(fovals(mat,100));
% grid on; box on;
% set(0,'DefaultFigureWindowStyle','docked')

% Compute condition number
% condestMat = condest(mat);
% 
% disp(['    Cond(mat)           ' num2str(condestMat,'%1.2e')]);
% disp(['---------------------------------------------------------']);

% eigtool(full(sysA.matA))

% -------------------------------------------------------------------------
% Compute iterative solution
% -------------------------------------------------------------------------

% solver = 'GMRES';
% switch solver
%     case 'CGNR'
%         [resRedVec, resPhyVec, errorVec] = solverCGNRredu_CG(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
%     case 'GMRES'
%         [resVec, ~, ~, ~, ~, X] = solverGMRES(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
% end

% b = sysA.rhsA;
% vecb = cat(1,b,-b);
rng("default")
vecb = rand(size(mat,2),1);
% vecb = zeros(size(mat,2),1)+1;

% [~,~,~,it,resVec] = gmres(mat,vecb,[],tol,maxit);
[~,~,~,it,resVec] = minres(mat,vecb,tol,maxit);

it_max = size(resVec,1);
bound = zeros(size(resVec,1),1);
for i=1:it_max    
    bound(i) = 2*((lambdaMax-lambdaMin)/(lambdaMax+lambdaMin))^(floor(i/2));
end

% solGMRES = X(:,end);
% writeField2D(dofm, mesh, solGMRES, 'output/solGMRES.pos', "solGMRES");
% writeField2D(dofm, mesh, solP, 'output/solInc.pos', "solInc");
% writeField2D(dofm, mesh, solGMRES+solP, 'output/solTot.pos', "solTot");

% get the last non zero residual
lastNonZero = find(resVec,1,'last');


figure
hold on
set(0,'DefaultFigureWindowStyle','docked')


% plot([0 maxit],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (direct)','linewidth', 1);
% plot([0 maxit],[errorProjL2 errorProjL2],'k:','DisplayName','Relative L2-error (projection)','linewidth', 1);
% semilogy(0:itout:maxit,resVec,'b-o','DisplayName','Relative residual','linewidth', 1,'markersize', 5);
semilogy(1:it_max,resVec,'b-o','DisplayName','Relative residual','linewidth', 1,'markersize', 5);
% semilogy(0:itout:maxit,errorVec,'k-o','DisplayName','Relative L2-error (iterative)','linewidth', 1,'markersize', 5);
semilogy(1:it_max,bound,'r--','DisplayName','Bound','linewidth', 1);
set(gca, 'YScale', 'log')
box on
grid on
xlim([0 lastNonZero+1]);
xlim([0 maxit]);
ylim auto;
% title(['CG - ' benchmark ' - ' solver ' - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)'], 'interpreter', 'latex', 'fontsize', 20)
xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
legend('Location', 'southwest', 'fontsize', 15)