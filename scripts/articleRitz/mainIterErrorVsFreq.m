%close all;
clear;

global k h

plotFlag = 1;

%% % CAVITY BENCHMARK
k = 3.01*sqrt(2)*pi;
h = 1/32;
% rangeFreq = 12.7:0.009:13.8;
rangeFreq = 12.7:0.1:13.8;
tol = 1e-10;
maxit = 2000;
L = 1;
degree = 2;
nbEigVec=11;


%% Without preconditioner
PREC = 0;
restart = 0;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);

%% With preconditioner
PREC = 1;
restart = 0;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 3;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 5;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);



%% % SCATTERING OPEN CAVITY DIRICHLET BENCHMARK
global LdomX LdomY LpmlX LpmlY
k = 23.676;
h = 1/20;
rangeFreq = 23.5:0.01:24.5;
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
restart = 0;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);

%% With preconditioner
PREC = 1;
restart = 0;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 3;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 5;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);



%% % SCATTERING OPEN CAVITY NEUMANN BENCHMARK
global LdomX LdomY LpmlX LpmlY
k = 23.591;
% k = 23.82;
% k = 24.275;
rangeFreq = 23.5:0.01:24.5;
h = 1/20;
tol = 1e-6;
maxit = 5000;
LdomX = 0.95;
LdomY = 0.5;
LpmlX = 0.2;
LpmlY = 0.2;
degree = 3;
nbEigVec=11;

%% Without preconditioner
PREC = 0;
restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);

%% With preconditioner
PREC = 1;
restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 3;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 5;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);


%%  %% MAIN FUNCTION %% %%


function run(benchmark,degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag)

    global k h


    mesh = setupBenchmark2D(benchmark);
    mesh = buildConnectivity2D(mesh);
    dofm = buildDofManager2D_CG(mesh, degree);

    if nbEigVec > 0
        nbit = zeros(length(rangeFreq),3);
        labels = ["k", "itG", "itD"];
    else
        nbit = zeros(length(rangeFreq),2);
        labels = ["k", "itG"];
    end

    nbit(:,1) = rangeFreq;

    switch benchmark
        case 'cavity'
            [eigvec,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,'firstEigvec',k);
            error = zeros(length(rangeFreq),2);
            error(:,1) = rangeFreq;
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

    for k = rangeFreq

        Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

        disp(['---------------------------------------------------------']);
        disp(['Method CG - Benchmark "' benchmark '"']);
        disp(['---------------------------------------------------------']);
        disp(['    k = ' num2str(k) ' - h = ' num2str(h) ' - Dlambda = ' num2str(Dlambda) ' - degree = ' num2str(degree) ' - PREC = ' num2str(PREC) ' - Restart = ' num2str(restart)]);
        disp(['---------------------------------------------------------']);

        [~, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);

        A = sysA.matA;
        M = sysA.matP;
        b = sysA.rhsA;
        AMinv = A/M;

        [Pdef,~,~] = computeDefOp(nbEigVec, eigvec, A);

        if k == rangeFreq(1)
            if maxit > size(A,2)
                maxit = size(A,2);
            end
    
            if restart == 0
                m = size(A,2);
            else
                m = restart;
                maxit = ceil(maxit/m);
            end
        end

        % No deflation

        % Compute GMRES solution
        disp(['|  GMRES - No deflation - PREC = ' num2str(PREC) ' - Restart = ' num2str(restart) '']);
        [uG, ~, ~, itG, ~] = gmres(AMinv,b,m,tol,maxit);
        itG = itG(2) + (itG(1)-1)*m;
        disp(['    converged in ' num2str(itG) ' iterations']);

        nbit(rangeFreq==k,2) = itG;

        if strcmp(benchmark, 'cavity')
            xG = M\uG;
            errorL2 = computeNormError2D_CG(mesh, dofm, xG);

            error(rangeFreq==k,2) = errorL2;
        end

        if nbEigVec > 0
            % Deflation

            % Compute GMRES solution
            disp(['|  GMRES - Deflation - PREC = ' num2str(PREC) ' - Restart = ' num2str(restart) '']);
            [~, ~, ~, itD, ~] = gmres(Pdef*AMinv,Pdef*b,m,tol,maxit);
            itD = itD(2) + (itD(1)-1)*m;
            disp(['    converged in ' num2str(itD) ' iterations']);

            nbit(rangeFreq==k,3) = itD;
        end



    end

    kmin = min(rangeFreq);
    kmax = max(rangeFreq);

    namefile = sprintf('output/iterVsFreq_%s_p%i_prec%i_range_%g-%g_def_%g_restart_%g.csv', benchmark, degree, PREC, kmin, kmax, nbEigVec, restart);
    writematrix([labels; nbit], namefile, 'Delimiter', 'comma');

    if strcmp(benchmark, 'cavity')
        namefile = sprintf('output/errorVsFreq_cavity_p%i_prec%i_range_%g-%g_def_%g_restart_%g.csv', degree, PREC, kmin, kmax, nbEigVec, restart);
        writematrix([["k", "errorL2"]; error], namefile, 'Delimiter', 'comma');
    end

    if plotFlag

        green = [0.4660 0.6740 0.1880];
        orange = [0.9290 0.6940 0.1250];

        figure;
        hold on;
        set(0,'DefaultFigureWindowStyle','docked');

        p1 = plot(nbit(:,1),nbit(:,2),'b-x','DisplayName',strcat('No deflation - PREC = ', num2str(PREC), ' - Restart = ', num2str(restart)),'linewidth', 2,'markersize', 10);
        p1.Color = green;
        if nbEigVec > 0
            p2 = plot(nbit(:,1),nbit(:,3),'r-o','DisplayName',strcat('Deflation (', num2str(nbEigVec), ' vec) - PREC = ', num2str(PREC), ' - Restart = ', num2str(restart)),'linewidth', 2,'markersize', 10);
            p2.Color = orange;
        end
        xlabel('k', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('Number of iterations', 'Interpreter', 'latex', 'FontSize', 15);
        title(['Nb of it vs freq - Benchmark "' benchmark '" - h = ' num2str(h) ' - degree = ' num2str(degree) ' - PREC = ' num2str(PREC) ' - Restart = ' num2str(restart)], 'Interpreter', 'latex', 'FontSize', 20);
        legend('Location', 'best', 'FontSize', 15);
    end

end