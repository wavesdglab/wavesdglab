%close all;
clear;
clear global;

global k h Options

Options.Basis = 'Jacobi'; % Jacobi, Lagrange
Options.Error = 'L2'; % L2, H1

plotFlag = 1;
saveSolFlag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Benchmark 'Cavity'
% Alternatives for PREC: 'none', 'CSL', 'ILU(CSL)', 'ILU(A)'

k = 3.01*sqrt(2)*pi;
h = 1/32;
tol = 1e-10;
degree = 2;

figure;
maxit = 300; itout = 5;
PREC = 'none'; nbEigVecList = [0, 1, 11]; restart = 0;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVecList,restart,plotFlag,saveSolFlag);

PREC = 'none'; nbEigVecList = [0, 1, 11]; restart = 25;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVecList,restart,plotFlag,saveSolFlag);

figure;
maxit = 100; itout = 1;
PREC = 'ILU(CSL)'; nbEigVecList = [0, 1, 11]; restart = 0;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVecList,restart,plotFlag,saveSolFlag);
PREC = 'ILU(CSL)'; nbEigVecList = [0, 1, 11]; restart = 25;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVecList,restart,plotFlag,saveSolFlag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Benchmark 'Scattering Open Cavity Neumann'
% Alternative with Dirichlet: 'scattering_openCavity_DIR' with 23.676

global LdomX LdomY LpmlX LpmlY WRITE_FIELD_ABSOLUTE
LdomX = 0.95; LdomY = 0.5; LpmlX = 0.2; LpmlY = 0.2;
WRITE_FIELD_ABSOLUTE = 1;

k = 23.591;
h = 1/20;
tol = 1e-6;
degree = 3;

figure;
maxit = 1450; itout = 10;
PREC = 'none'; nbEigVec=[0, 1, 10]; restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
PREC = 'none'; nbEigVec=[0, 1, 10]; restart = 25;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);

figure;
maxit = 1000; itout = 10;
PREC = 'ILU(CSL)'; nbEigVec=[0, 1, 10]; restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
PREC = 'ILU(CSL)'; nbEigVec=[0, 1, 10]; restart = 250;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,PREC,tol,maxit,itout,nbEigVecList,restart,plotFlag,saveSolFlag)

global k h

mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

disp(['---------------------------------------------------------']);
disp(['Method CG - Benchmark "' benchmark '"']);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    PREC                ' PREC]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

switch PREC
    case 'none'
        prec = 0;
    case {'CSL','ILU(CSL)','ILU(A)'}
        prec = 1;
    otherwise
        error('Error. \n%s is not a valid preconditioner', PREC);
end

[~, sysA] = computeSolNum2D_CG(mesh, dofm, prec);

A = sysA.matA;
P = sysA.matP;
b = sysA.rhsA;

invPrec = @(x) P\x;
switch PREC
    case 'ILU(CSL)'
        [L,U] = ilu(P);
        invPrec = @(x) U\(L\x);
    case 'ILU(A)'
        [L,U] = ilu(A);
        invPrec = @(x) U\(L\x);
end

maxit = min(maxit, size(A,2));

if restart == 0
    m = size(A,2);
else
    m = restart;
    maxit = ceil(maxit/m);
end

for iterDeflation = 1:size(nbEigVecList(:),1)
    nbEigVec = nbEigVecList(iterDeflation);

    disp(['|  GMRES - PREC = ' PREC ' - Deflation = ' num2str(nbEigVec) ' - Restart = ' num2str(restart) '']);

    % Compute deflation matrices
    if nbEigVec > 0
        switch benchmark
            case 'cavity'
                [eigvec,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,'closestEigvec',k);
            case 'scattering_openCavity_NEU'
                [eigvec,nbEigVec] = computeProjEigVec_openCavity_NEU(mesh, dofm, nbEigVec, k);
            case 'scattering_openCavity_DIR'
                [eigvec,nbEigVec] = computeProjEigVec_openCavity_DIR(mesh, dofm, nbEigVec, k);
            otherwise
                error('Error. \n%s is not a valid benchmark for deflation', benchmark);
        end
        [PdefA,Pdef,Qdef,Q] = computeDefOp(nbEigVec, eigvec, A);
    else
        PdefA = @(x) A*x;
        Pdef = 1;
        Qdef = 1;
        Q = 0;
    end

    % Compute GMRES solution
    [uD, ~, ~, it, vecRes] = gmres(@(x) PdefA(invPrec(x)), Pdef*b, m, tol, maxit);
    disp(['    converged in ' num2str(size(vecRes,1)-1) ' iterations']);

    vecRes = vecRes(:)./vecRes(1);
    vecRes = vecRes(1:itout:end);
    vecIter = 0:itout:itout*size(vecRes,1)-1;
    xD = Q*b + Qdef*(invPrec(uD));

    if saveSolFlag
        namefile = sprintf('output/numSolD_%s_p%i_k=%g_prec=%s_def=%g_restart=%g.pos', benchmark, degree, k, PREC, nbEigVec, restart);
        namesol = strcat('x_k=', num2str(k), '_PREC=', PREC, '_def=', num2str(nbEigVec), '_restart=', num2str(restart));
        writeField2D(dofm, mesh, xD, namefile, namesol);
    end

    namefile = sprintf('output/historyGMRES_%s_p%i_k=%g_prec=%s_def=%g_restart=%g.csv', benchmark, degree, k, PREC, nbEigVec, restart);
    writematrix([["it", "rrG"]; [vecIter' vecRes]], namefile, 'Delimiter', 'comma');

    if plotFlag
        hold on;
        semilogy(vecIter, vecRes, 'DisplayName',['PREC = ' PREC ' - nDef = ' num2str(nbEigVec) ' - nRestart = ' num2str(restart)]);
        set(gca, 'YScale', 'log')
        box on;
        grid on;
        ylim auto;
        xlabel('Iteration number');
        ylabel('Relative residual');
        title(['Benchmark "' benchmark '" - k = ' num2str(k) ' - h = ' num2str(h) ' - degree = ' num2str(degree) ' - Dlambda = ' num2str(Dlambda)]);
        legend('Location', 'southwest');
        ylim([0 maxit]);
        ylim([tol 1]);
    end

end

% These variables must be cleared if 'cavity' is run after 'scattering'.
clear global LdomX LdomY LpmlX LpmlY;

end