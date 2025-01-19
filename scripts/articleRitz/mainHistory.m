%close all;
clear;

global k h

plotFlag = 1;
saveSolFlag = 0;

%% % CAVITY BENCHMARK
k = 3.01*sqrt(2)*pi;
h = 1/32;
tol = 1e-10;
maxit = 2000;
L = 1;
degree = 2;
nbEigVec=1;


%% Without preconditioner
PREC = 0;
itout =4;
restart = 0;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 10;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 25;
run('cavity',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);

%% With preconditioner
PREC = 1;
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
PREC = 0;
itout = 10;
restart = 0;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 10;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 25;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);

%% With preconditioner
PREC = 1;
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
PREC = 0;
itout = 10;
restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 10;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);
restart = 25;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag);

%% With preconditioner
PREC = 1;
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
    disp(['    PREC                ' num2str(PREC)]);
    disp(['    Dlambda             ' num2str(Dlambda)]);
    disp(['---------------------------------------------------------']);

    [~, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);

    A = sysA.matA;
    M = sysA.matP;
    b = sysA.rhsA;
    AMinv = A/M;

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
    disp(['|  GMRES - No deflation - PREC = ' num2str(PREC) ' - Restart = ' num2str(restart) '']);
    [uG, ~, ~, itG, rrG] = gmres(AMinv,b,m,tol,maxit);
    itG = itG(2) + (itG(1)-1)*m;
    rrG = rrG(:)./rrG(1);
    rrG = rrG(1:itout:end);
    iterG = 0:itout:itout*size(rrG,1)-1;
    rrG = [iterG' rrG];
    xG = M\uG;
    disp(['    converged in ' num2str(itG) ' iterations']);

    if saveSolFlag
        namefile = sprintf('output/numSolG_%s_p%i_prec%i_k_%g_def_%g_restart_%g.pos', benchmark, degree, PREC, k, nbEigVec, restart);
        writeField2D(dofm, mesh, xG, namefile, "xG");
    end

    labels = ["it", "rrG"];
    results = rrG;

    if nbEigVec > 0
        % Deflation

        % Compute GMRES solution
        disp(['|  GMRES - Deflation - PREC = ' num2str(PREC) ' - Restart = ' num2str(restart) '']);
        [uD, ~, ~, itD, rrD] = gmres(Pdef*AMinv,Pdef*b,m,tol,maxit);
        itD = itD(2) + (itD(1)-1)*m;
        rrD = rrD(:)./rrD(1);
        rrD = rrD(1:itout:end);
        iterD = 0:itout:itout*size(rrD,1)-1;
        rrD = [iterD' rrD];
        xD = Q*b + Qdef*(M\uD);
        disp(['    converged in ' num2str(itD) ' iterations']);

        if saveSolFlag
            namefile = sprintf('output/numSolD_%s_p%i_prec%i_k_%g_def_%g_restart_%g.pos', benchmark, degree, PREC, k, nbEigVec, restart);
            writeField2D(dofm, mesh, xD, namefile, "xD");
        end

        labels = ["it", "rrG", "rrD"];
        % Add zeros to rrD to match the size of rrG
        rrD = [rrD; zeros(size(rrG,1)-size(rrD,1),2)];
        results = [rrG rrD(:,2)];
    end

    namefile = sprintf('output/historyGMRES_%s_p%i_prec%i_k_%g_def_%g_restart_%g.csv', benchmark, degree, PREC, k, nbEigVec, restart);
    writematrix([labels; results], namefile, 'Delimiter', 'comma');

    if plotFlag

        green = [0.4660 0.6740 0.1880];
        orange = [0.9290 0.6940 0.1250];

        figure;
        hold on
        set(0,'DefaultFigureWindowStyle','docked')

        p1 = semilogy(rrG(:,1), rrG(:,2), 'b-x','DisplayName',strcat('No deflation - PREC = ', num2str(PREC), ' - Restart = ', num2str(restart)),'linewidth', 2,'markersize', 10);
        p1.Color = green;
        if nbEigVec > 0
            p2 = semilogy(rrD(:,1), rrD(:,2), 'r-o','DisplayName',strcat('Deflation (', num2str(nbEigVec), ' vec) - PREC = ', num2str(PREC), ' - Restart = ', num2str(restart)),'linewidth', 2,'markersize', 10);
            p2.Color = orange;
        end
        set(gca, 'YScale', 'log')
        box on;
        grid on;
        ylim auto;
        xlabel('Iteration number', 'interpreter', 'latex', 'fontsize', 15);
        ylabel('Relative residual', 'interpreter', 'latex', 'fontsize', 15);
        title(['GMRES - "' benchmark '" - k = ' num2str(k) ' - h = ' num2str(h) ' - degree = ' num2str(degree) ' - PREC = ' num2str(PREC) ' - Dlambda = ' num2str(Dlambda) ' - Restart = ' num2str(restart)], 'interpreter', 'latex', 'fontsize', 20);
        legend('Location', 'southwest', 'fontsize', 15);
    end

end