%close all;
clear;
clear global;

N=5;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);

global k h

plotFlag = 1;
saveSolFlag = 1;

%% % CAVITY BENCHMARK
k = 3.01*sqrt(2)*pi;
h = 1/32;
tol = 1e-10;
maxit = 2000;
L = 1;
degree = 2;
nbEigVec=1;


% %% Without preconditioner
PREC = 'none';
itout =4;
restart = 0;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 10;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 25;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);

%% With preconditioner: CSL
PREC = 'CSL';
itout = 1;
restart = 0;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 3;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 5;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);

%% With preconditioner: CSL (ILU)
PREC = 'CSL(ILU)';
itout = 1;
restart = 0;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 3;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 5;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);



%% % SCATTERING OPEN CAVITY DIRICHLET BENCHMARK
global LdomX LdomY LpmlX LpmlY
k = 23.676;
h = 1/20;
tol = 1e-6;
maxit = 5000;
LdomX = 0.95;
LdomY = 0.5;
LpmlX = 0.2;
LpmlY = 0.2;
degree = 3;
nbEigVec = 1;

%% Without preconditioner
PREC = 'none';
itout = 10;
restart = 0;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 10;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 25;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);

%% With preconditioner: CSL
PREC = 'CSL';
itout = 1;
restart = 0;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 3;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 5;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);

%% With preconditioner: CSL (ILU)
PREC = 'CSL(ILU)';
itout = 1;
restart = 0;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 3;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 5;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);



%% % SCATTERING OPEN CAVITY NEUMANN BENCHMARK
global LdomX LdomY LpmlX LpmlY
k = 23.591;
% k = 23.82;
% k = 24.275;
h = 1/20;
tol = 1e-6;
maxit = 5000;
LdomX = 0.95;
LdomY = 0.5;
LpmlX = 0.2;
LpmlY = 0.2;
degree = 3;
nbEigVec=1;

%% Without preconditioner
PREC = 'none';
itout = 10;
restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 10;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 25;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);

%% With preconditioner: CSL
PREC = 'CSL';
itout = 1;
restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 3;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 5;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);

%% With preconditioner: CSL (ILU)
PREC = 'CSL(ILU)';
itout = 1;
restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 3;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 5;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);

%%  %% MAIN FUNCTION %% %%


function run(benchmark,degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag)

    global k h
    
    if strcmp(benchmark, 'cavity')
        clear global LdomX LdomY LpmlX LpmlY
    end

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
        case 'CSL'
            prec = 1;
        case 'CSL(ILU)'
            prec = 1;
        otherwise
            error('Unknown preconditioner');
    end

    [~, sysA] = computeSolNum2D_CG(mesh, dofm, prec);

    A = sysA.matA;
    M = sysA.matP;
    b = sysA.rhsA;

    invPrec = @(x) M\x;

    if strcmp(PREC, 'CSL(ILU)')
        [L,U] = ilu(M);
        invPrec = @(x) U\(L\x);
    end


    switch benchmark
        case 'cavity'
            [eigvec,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,'closestEigvec',k);
        case 'scattering_openCavity_NEU'
            [eigvec,nbEigVec] = computeProjEigVec_openCavity_NEU(mesh, dofm, nbEigVec, k);
            global WRITE_FIELD_ABSOLUTE
            WRITE_FIELD_ABSOLUTE = 1;
        case 'scattering_openCavity_DIR'
            [eigvec,nbEigVec] = computeProjEigVec_openCavity_DIR(mesh, dofm, nbEigVec, k);
            global WRITE_FIELD_ABSOLUTE
            WRITE_FIELD_ABSOLUTE = 1;
        otherwise
            error('Unknown benchmark');
    end

    [Pdef,Qdef,Q] = computeDefOp(nbEigVec, eigvec, A);

    if maxit > size(A,2)
        maxit = size(A,2);
    end

    if restart == 0
        m = size(A,2);
    else
        m = restart;
        maxit = ceil(maxit/m);
    end

    % No deflation

    % Compute GMRES solution
    disp(['|  GMRES - No deflation - PREC = ' PREC ' - Restart = ' num2str(restart) '']);
    [uG, ~, ~, itG, rrG] = gmres(@(x) A*(invPrec(x)),b,m,tol,maxit);
    itG = itG(2) + (itG(1)-1)*m;
    rrG = rrG(:)./rrG(1);
    rrG = rrG(1:itout:end);
    iterG = 0:itout:itout*size(rrG,1)-1;
    rrG = [iterG' rrG];
    xG = invPrec(uG);
    disp(['    converged in ' num2str(itG) ' iterations']);

    if saveSolFlag
        namefile = sprintf('output/numSolG_%s_p%i_k=%g_prec=%s_def=%g_restart=%g.pos', benchmark, degree, k, PREC, nbEigVec, restart);
        namesol = strcat('x_k=', num2str(k), '_PREC=', PREC, '_restart=', num2str(restart));
        writeField2D(dofm, mesh, xG, namefile, namesol);
    end

    labels = ["it", "rrG"];
    results = rrG;

    if nbEigVec > 0
        % Deflation

        % Compute GMRES solution
        disp(['|  GMRES - Deflation - PREC = ' PREC ' - Restart = ' num2str(restart) '']);
        [uD, ~, ~, itD, rrD] = gmres(@(x) Pdef*(A*(invPrec(x))),Pdef*b,m,tol,maxit);
        itD = itD(2) + (itD(1)-1)*m;
        rrD = rrD(:)./rrD(1);
        rrD = rrD(1:itout:end);
        iterD = 0:itout:itout*size(rrD,1)-1;
        rrD = [iterD' rrD];
        xD = Q*b + Qdef*(invPrec(uD));
        disp(['    converged in ' num2str(itD) ' iterations']);

        if saveSolFlag
            namefile = sprintf('output/numSolD_%s_p%i_k=%g_prec=%s_def=%g_restart=%g.pos', benchmark, degree, k, PREC, nbEigVec, restart);
            namesol = strcat('x_k=', num2str(k), '_PREC=', PREC, '_def=', num2str(nbEigVec), '_restart=', num2str(restart));
            writeField2D(dofm, mesh, xD, namefile, namesol);
        end

        labels = ["it", "rrG", "rrD"];
        rrD = [rrD; zeros(size(rrG,1)-size(rrD,1),2)];
        results = [rrG rrD(:,2)];
    end

    namefile = sprintf('output/historyGMRES_%s_p%i_k_%g_prec=%s_def_%g_restart_%g.csv', benchmark, degree, k, PREC, nbEigVec, restart);
    writematrix([labels; results], namefile, 'Delimiter', 'comma');

    if plotFlag

        green = [0.4660 0.6740 0.1880];
        orange = [0.9290 0.6940 0.1250];

        figure;
        hold on
        set(0,'DefaultFigureWindowStyle','docked')

        p1 = semilogy(rrG(:,1), rrG(:,2), 'b-x','DisplayName',strcat('No deflation - PREC = ', PREC, ' - Restart = ', num2str(restart)),'linewidth', 2,'markersize', 10);
        p1.Color = green;
        if nbEigVec > 0
            p2 = semilogy(rrD(:,1), rrD(:,2), 'r-o','DisplayName',strcat('Deflation (', num2str(nbEigVec), ' vec) - PREC = ', PREC, ' - Restart = ', num2str(restart)),'linewidth', 2,'markersize', 10);
            p2.Color = orange;
        end
        set(gca, 'YScale', 'log')
        box on;
        grid on;
        ylim auto;
        xlabel('Iteration number', 'interpreter', 'latex', 'fontsize', 15);
        ylabel('Relative residual', 'interpreter', 'latex', 'fontsize', 15);
        title(['GMRES - "' benchmark '" - k = ' num2str(k) ' - h = ' num2str(h) ' - degree = ' num2str(degree) ' - PREC = ', PREC, ' - Dlambda = ' num2str(Dlambda) ' - Restart = ' num2str(restart)], 'interpreter', 'latex', 'fontsize', 20);
        legend('Location', 'southwest', 'fontsize', 15);
    end

end