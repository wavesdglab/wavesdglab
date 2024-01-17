% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Authors: Axel Modave, Pierre Marchand, Timothée Raynaud

clear all;
%close all;

global k L L_PML R_disk;
L_PML = 0;
R_disk = 0;
N=10;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);


% Setup benchmark and parameters
benchmark = 'cavity';
switch benchmark
    case 'open'
        k = 15*pi;
        h = 1/16;
        tol = 1e-10; maxit = 1000; itout = 50;
    case 'cavity'
        k = 0.1*sqrt(2)*pi;
        h = 1/16;
        tol = 1e-10; maxit = 2000; itout =100;
        L = 1;
    case 'scatteringPML'
        k = 15;
        h = 0.02;
        tol = 1e-10; maxit = 2000; itout = 1000;
        L = 1.1;
        R_disk = 1;
        L_PML = 0.2;
    case 'waveguide'
        k = 6*pi;
        h = 1/8;
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 1; % P1
PREC = 1; % for preconditioner

% Build mesh and DOF manager
mesh = benchmark2D(benchmark,h);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree); % espace fonctionnel discret

% tab_k = [2.1 2.01 2.0]*sqrt(2)*pi;
tab_tol = [1e-1 1e-2 1e-3 1e-4 1e-5];
k = 0.05;
k_max = 28;
imax = k_max/0.05;
cond_array = zeros(imax+1,4);
it_A = zeros(imax+1,6);
err_A = zeros(imax+1,6);
% it_B_left = zeros(imax+1,6);
% err_B_left = zeros(imax+1,6);
% it_C_left = zeros(imax+1,6);
% err_C_left = zeros(imax+1,6);

i = 1; % < imax + 1
while(k<k_max & i<imax+1)
% for i=1:length(tab_k)
    % k = tab_k(i);

    Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1); % nb de points par longueur d'onde

    % -------------------------------------------------------------------------
    % Compute solution and error
    % -------------------------------------------------------------------------

    disp(['---------------------------------------------------------']);
    disp(['Method CG - Benchmark "' benchmark '"']);
    disp(['---------------------------------------------------------']);
    % disp(['    k/(sqrt(2)*pi)      ' num2str(k/(sqrt(2)*pi))]);
    disp(['    k                   ' num2str(k)]);
    disp(['    k^2/(pi^2)          ' num2str(k^2/(pi^2))]);
    disp(['    h                   ' num2str(h)]);
    disp(['    degree              ' num2str(degree)]);
    disp(['    Dlambda             ' num2str(Dlambda)]);
    disp(['---------------------------------------------------------']);

    [sysA, ~] = computeSolNum2D_CG(mesh, dofm, PREC);
    % [sysB, ~] = computeSolNum2D_CG(mesh, dofm, PREC);
    % [sysC, ~] = computeSolNum2D_CG(mesh, dofm, PREC);

    % errorL2 = computeNormError2D_CG(mesh, dofm, solA);

    % solP = computeSolProjL2_2D_CG(mesh, dofm);
    % errorProjL2 = computeNormError2D_CG(mesh, dofm, solP);
 
    % disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
    % disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
    % disp(['---------------------------------------------------------']);

    matA = sysA.matA; % A
    % matprecB = sysB.matP\sysB.matA; % M^-1 A
    % matprecC = sysC.matP\sysC.matA; % SL^-1 A

    condmatA = condest(matA);
    % condmatB = condest(matprecB);
    % condmatC = condest(matprecC);
    
    cond_array(i,1) = k;
    cond_array(i,2) = condmatA;
    % cond_array(i,3) = condmatB;
    % cond_array(i,4) = condmatC;

    % compute solverGMRES for different tolerance
    disp(['---------------------------------------------------------']);
    disp(['NO PREC ']);
    [~, ~, it, flag, err, x] = solverGMRES(mesh, dofm, sysA, tab_tol, maxit, itout, @computeNormError2D_CG);

    x1 = x(:,1);
    x2 = x(:,2);
    x3 = x(:,3);
    x4 = x(:,4);
    x5 = x(:,5);

    it_A(i,1) = k;
    it_A(i,2:end) = it;
    err_A(i,1) = k;
    err_A(i,2:end) = err;

    % disp(['---------------------------------------------------------']);
    % disp(['PREC LEFT BY M^-1 ']);
    % [~, ~, it, flag, err, x] = solverGMRES_left(mesh, dofm, sysB, tab_tol, maxit, itout, @computeNormError2D_CG);

    % x1 = x(:,1);
    % x2 = x(:,2);
    % x3 = x(:,3);
    % x4 = x(:,4);
    % x5 = x(:,5);

    % it_B_left(i,1) = k;
    % it_B_left(i,2:end) = it;
    % err_B_left(i,1) = k;
    % err_B_left(i,2:end) = err;


    % disp(['---------------------------------------------------------']);
    % disp(['PREC LEFT BY SHIFTED LAPLACIAN ']);
    % [~, ~, it, flag, err, x] = solverGMRES_left(mesh, dofm, sysC, tab_tol, maxit, itout, @computeNormError2D_CG);

    % x1 = x(:,1);
    % x2 = x(:,2);
    % x3 = x(:,3);
    % x4 = x(:,4);
    % x5 = x(:,5);

    % it_C_left(i,1) = k;
    % it_C_left(i,2:end) = it;
    % err_C_left(i,1) = k;
    % err_C_left(i,2:end) = err;
    
    % k = k + 0.1*sqrt(2)*pi;
    k = k + 0.05;
    i = i + 1;

end

csvwrite('output/it_A.csv', it_A);
csvwrite('output/err_A.csv', err_A);
% csvwrite('output/it_B_left.csv', it_B_left);
% csvwrite('output/err_B_left.csv', err_B_left);
% csvwrite('output/it_C_left.csv', it_C_left);
% csvwrite('output/err_C_left.csv', err_C_left);