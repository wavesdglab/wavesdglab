%close all;
clear;

N=10;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);

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
nbEigVec=11;
itout = 1;
numTests = 10;

restart = 0;
run('cavity',degree,tol,maxit,itout,nbEigVec,restart,numTests,plotFlag,saveSolFlag);
restart = 3;
run('cavity',degree,tol,maxit,itout,nbEigVec,restart,numTests,plotFlag,saveSolFlag);
restart = 5;
run('cavity',degree,tol,maxit,itout,nbEigVec,restart,numTests,plotFlag,saveSolFlag);



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
itout = 1;
numTests = 10;

restart = 0;
run('scattering_openCavity_DIR',degree,tol,maxit,itout,nbEigVec,restart,numTests,plotFlag,saveSolFlag);
restart = 3;
run('scattering_openCavity_DIR',degree,tol,maxit,itout,nbEigVec,restart,numTests,plotFlag,saveSolFlag);
restart = 5;
run('scattering_openCavity_DIR',degree,tol,maxit,itout,nbEigVec,restart,numTests,plotFlag,saveSolFlag);

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
itout = 1;
numTests = 10;

restart = 0;
run('scattering_openCavity_NEU',degree,tol,maxit,itout,nbEigVec,restart,numTests,plotFlag,saveSolFlag);
restart = 3;
run('scattering_openCavity_NEU',degree,tol,maxit,itout,nbEigVec,restart,numTests,plotFlag,saveSolFlag);
restart = 5;
run('scattering_openCavity_NEU',degree,tol,maxit,itout,nbEigVec,restart,numTests,plotFlag,saveSolFlag);

%%  %% MAIN FUNCTION %% %%


function run(benchmark,degree,tol,maxit,itout,nbEigVec,restart,numTests,plotFlag,saveSolFlag)

    global k h

    PREC = 1;
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
    disp(['    Dlambda             ' num2str(Dlambda)]);
    disp(['---------------------------------------------------------']);

    [~, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);

    A = sysA.matA;
    M = sysA.matP;
    b = sysA.rhsA;

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

    labelsGMRES = ["No def - No ILU", "No def - ILU(0)", "No def - ILUC", "No def - ILUTP"];
    labelsILU = ["ILU(0)", "ILUC", "ILUTP"];
    timeGMRES = zeros(numTests+3, length(labelsGMRES));
    timeLU = zeros(numTests, length(labelsILU));
    errL2 = zeros(1, length(labelsGMRES));
    it = zeros(1, length(labelsGMRES));

    for l = 1:numTests
        tic;
        [uG, ~, ~, itG, ~] = gmres(@(x) A*(M\x),b,m,tol,maxit);
        timeGMRES(l,1) = toc;
        if l==1
            xRef = M\uG;
            errL2(1) = 0;
            itG = itG(2) + (itG(1)-1)*m;
            it(1) = itG;
        end
    end

    timeGMRES(numTests+1,1) = mean(timeGMRES(1:numTests,1));
    timeGMRES(numTests+2,1) = std(timeGMRES(1:numTests,1));
    timeGMRES(numTests+3,1) = timeGMRES(numTests+1,1)/it(1);

    disp(['|  No deflation - No ILU - iterations   = ' num2str(it(1))]);
    disp(['|                          Mean time    = ' num2str(timeGMRES(numTests+1,1)) 's']);
    disp(['|                          Std time     = ' num2str(timeGMRES(numTests+2,1)) 's']);
    disp(['|                          Mean time/it = ' num2str(timeGMRES(numTests+3,1)) 's']);

    for l = 1:numTests
        tic;
        [L,U] = ilu(M);
        timeLU(l,1) = toc;
        tic;
        [uG, ~, ~, itG, ~] = gmres(@(x) A*(U\(L\x)),b,m,tol,maxit);
        timeGMRES(l,2) = toc;
        if l==1
            x = U\(L\uG);
            errL2(2) = norm(xRef-x);
            itG = itG(2) + (itG(1)-1)*m;
            it(2) = itG;
        end
    end

    timeLU(numTests,1) = mean(timeLU(1:numTests,1));
    timeLU(numTests+1,1) = std(timeLU(1:numTests,1));
    timeGMRES(numTests+1,2) = mean(timeGMRES(1:numTests,2));
    timeGMRES(numTests+2,2) = std(timeGMRES(1:numTests,2));
    timeGMRES(numTests+3,2) = timeGMRES(numTests+1,2)/it(2);

    disp(['|  No deflation - ILU(0) - iterations      = ' num2str(it(2))]);
    disp(['|                          Mean time       = ' num2str(timeGMRES(numTests+1,2)) 's']);
    disp(['|                          Std time        = ' num2str(timeGMRES(numTests+2,2)) 's']);
    disp(['|                          Mean time/it    = ' num2str(timeGMRES(numTests+3,2)) 's']);
    disp(['|                          Mean LU(0) time = ' num2str(timeLU(numTests,1)) 's']);
    disp(['|                          Std LU(0) time  = ' num2str(timeLU(numTests+1,1)) 's']);
    disp(['|                          ||xRef - x||    = ' num2str(errL2(2))]);

    for l = 1:numTests
        tic;
        [L,U] = ilu(M,struct('type','crout','droptol',1e-3));
        timeLU(l,2) = toc;
        tic;
        [uG, ~, ~, itG, ~] = gmres(@(x) A*(U\(L\x)),b,m,tol,maxit);
        timeGMRES(l,3) = toc;
        if l==1
            x = U\(L\uG);
            errL2(3) = norm(xRef-x);
            itG = itG(2) + (itG(1)-1)*m;
            it(3) = itG;
        end
    end

    timeLU(numTests,2) = mean(timeLU(1:numTests,2));
    timeLU(numTests+1,2) = std(timeLU(1:numTests,2));
    timeGMRES(numTests+1,3) = mean(timeGMRES(1:numTests,3));
    timeGMRES(numTests+2,3) = std(timeGMRES(1:numTests,3));
    timeGMRES(numTests+3,3) = timeGMRES(numTests+1,3)/it(3);

    disp(['|  No deflation - ILUC - iterations      = ' num2str(it(3))]);
    disp(['|                        Mean time       = ' num2str(timeGMRES(numTests+1,3)) 's']);
    disp(['|                        Std time        = ' num2str(timeGMRES(numTests+2,3)) 's']);
    disp(['|                        Mean time/it    = ' num2str(timeGMRES(numTests+3,3)) 's']);
    disp(['|                        Mean ILUC time  = ' num2str(timeLU(numTests,2)) 's']);
    disp(['|                        Std ILUC time   = ' num2str(timeLU(numTests+1,2)) 's']);
    disp(['|                        ||xRef - x||    = ' num2str(errL2(3))]);

    for l = 1:numTests
        tic;
        [L,U] = ilu(M,struct('type','ilutp','droptol',1e-3));
        timeLU(l,3) = toc;
        tic;
        [uG, ~, ~, itG, ~] = gmres(@(x) A*(U\(L\x)),b,m,tol,maxit);
        timeGMRES(l,4) = toc;
        if l==1
            x = U\(L\uG);
            errL2(4) = norm(xRef-x);
            itG = itG(2) + (itG(1)-1)*m;
        end
    end

    timeLU(numTests,3) = mean(timeLU(1:numTests,3));
    timeLU(numTests+1,3) = std(timeLU(1:numTests,3));
    timeGMRES(numTests+1,4) = mean(timeGMRES(1:numTests,4));
    timeGMRES(numTests+2,4) = std(timeGMRES(1:numTests,4));
    timeGMRES(numTests+3,4) = timeGMRES(numTests+1,4)/it(4);

    disp(['|  No deflation - ILUTP - iterations       = ' num2str(it(4))]);
    disp(['|                         Mean time        = ' num2str(timeGMRES(numTests+1,4)) 's']);
    disp(['|                         Std time         = ' num2str(timeGMRES(numTests+2,4)) 's']);
    disp(['|                         Mean time/it     = ' num2str(timeGMRES(numTests+3,4)) 's']);
    disp(['|                         Mean ILUTP time  = ' num2str(timeLU(numTests,3)) 's']);
    disp(['|                         Std ILUTP time   = ' num2str(timeLU(numTests+1,3)) 's']);
    disp(['|                         ||xRef - x||     = ' num2str(errL2(4))]);


    if nbEigVec > 0
        timeGMRES = [timeGMRES zeros(numTests+3,4)];
        labelsGMRES = ["No def - No ILU", "No def - ILU(0)", "No def - ILUC", "No def - ILUTP", "Def - No ILU", "Def - ILU(0)", "Def - ILUC", "Def - ILUTP"];
        errL2 = [errL2 zeros(1,4)];
        it = [it zeros(1,4)];
        % Deflation

        for l = 1:numTests
            tic;
            [uD, ~, ~, itD, ~] = gmres(@(x) Pdef*(A*(M\x)),Pdef*b,m,tol,maxit);
            timeGMRES(l,5) = toc;
            if l==1
                xD = Q*b + Qdef*(M\uD);
                errL2(5) = norm(xRef-xD);
                itD = itD(2) + (itD(1)-1)*m;
                it(5) = itD;
            end
        end

        timeGMRES(numTests+1,5) = mean(timeGMRES(1:numTests,5));
        timeGMRES(numTests+2,5) = std(timeGMRES(1:numTests,5));
        timeGMRES(numTests+3,5) = timeGMRES(numTests+1,5)/it(5);
        
        disp(['|  Deflation - No ILU - iterations   = ' num2str(it(5))]);
        disp(['|                       Mean time    = ' num2str(timeGMRES(numTests+1,5)) 's']);
        disp(['|                       Std time     = ' num2str(timeGMRES(numTests+2,5)) 's']);
        disp(['|                       Mean time/it = ' num2str(timeGMRES(numTests+3,5)) 's']);
        disp(['|                       ||xRef - x|| = ' num2str(errL2(5))]);

        [L,U] = ilu(M);
        for l = 1:numTests
            [uD, ~, ~, itD, ~] = gmres(@(x) Pdef*(A*(U\(L\x))),Pdef*b,m,tol,maxit);
            timeGMRES(l,6) = toc;
            if l==1
                x = Q*b + Qdef*(U\(L\uD));
                errL2(6) = norm(xRef-x);
                itD = itD(2) + (itD(1)-1)*m;
                it(6) = itD;
            end
        end

        timeGMRES(numTests+1,6) = mean(timeGMRES(1:numTests,6));
        timeGMRES(numTests+2,6) = std(timeGMRES(1:numTests,6));
        timeGMRES(numTests+3,6) = timeGMRES(numTests+1,6)/it(6);

        disp(['|  Deflation - ILU(0) - iterations      = ' num2str(it(6))]);
        disp(['|                       Mean time       = ' num2str(timeGMRES(numTests+1,6)) 's']);
        disp(['|                       Std time        = ' num2str(timeGMRES(numTests+2,6)) 's']);
        disp(['|                       Mean time/it    = ' num2str(timeGMRES(numTests+3,6)) 's']);
        disp(['|                       ||xRef - x||    = ' num2str(errL2(6))]);

        [L,U] = ilu(M,struct('type','crout','droptol',1e-3));
        for l = 1:numTests
            [uD, ~, ~, itD, ~] = gmres(@(x) Pdef*(A*(U\(L\x))),Pdef*b,m,tol,maxit);
            timeGMRES(l,7) = toc;
            if l==1
                x = Q*b + Qdef*(U\(L\uD));
                errL2(7) = norm(xRef-x);
                itD = itD(2) + (itD(1)-1)*m;
                it(7) = itD;
            end
        end

        timeGMRES(numTests+1,7) = mean(timeGMRES(1:numTests,7));
        timeGMRES(numTests+2,7) = std(timeGMRES(1:numTests,7));
        timeGMRES(numTests+3,7) = timeGMRES(numTests+1,7)/it(7);

        disp(['|  Deflation - ILUC - iterations      = ' num2str(it(7))]);
        disp(['|                     Mean time       = ' num2str(timeGMRES(numTests+1,7)) 's']);
        disp(['|                     Std time        = ' num2str(timeGMRES(numTests+2,7)) 's']);
        disp(['|                     Mean time/it    = ' num2str(timeGMRES(numTests+3,7)) 's']);
        disp(['|                     ||xRef - x||    = ' num2str(errL2(7))]);

        [L,U] = ilu(M,struct('type','ilutp','droptol',1e-3));
        for l = 1:numTests
            [uD, ~, ~, itD, ~] = gmres(@(x) Pdef*(A*(U\(L\x))),Pdef*b,m,tol,maxit);
            timeGMRES(l,8) = toc;
            if l==1
                x = Q*b + Qdef*(U\(L\uD));
                errL2(8) = norm(xRef-x);
                itD = itD(2) + (itD(1)-1)*m;
            end
        end

        timeGMRES(numTests+1,8) = mean(timeGMRES(1:numTests,8));
        timeGMRES(numTests+2,8) = std(timeGMRES(1:numTests,8));
        timeGMRES(numTests+3,8) = timeGMRES(numTests+1,8)/it(8);
     
        disp(['|  Deflation - ILUTP - iterations       = ' num2str(it(8))]);
        disp(['|                      Mean time        = ' num2str(timeGMRES(numTests+1,8)) 's']);
        disp(['|                      Std time         = ' num2str(timeGMRES(numTests+2,8)) 's']);
        disp(['|                      Mean time/it     = ' num2str(timeGMRES(numTests+3,8)) 's']);
        disp(['|                      ||xRef - x||     = ' num2str(errL2(8))]);

    end


    % namefile = sprintf('output/timeGMRES_%s_p%i_k_%g_def_%g_restart_%g.mat', benchmark, degree, k, nbEigVec, restart);
    % writematrix([labelsGMRES; timeGMRES], namefile, 'Delimiter', 'comma');

    % namefile = sprintf('output/timeLU_%s_p%i_k_%g_def_%g_restart_%g.mat', benchmark, degree, k, nbEigVec, restart);
    % writematrix([labelsILU; timeLU], namefile, 'Delimiter', 'comma');

    % namefile = sprintf('output/errL2_%s_p%i_k_%g_def_%g_restart_%g.mat', benchmark, degree, k, nbEigVec, restart);
    % writematrix(errL2, namefile, 'Delimiter', 'comma');

    % namefile = sprintf('output/it_%s_p%i_k_%g_def_%g_restart_%g.mat', benchmark, degree, k, nbEigVec, restart);
    % writematrix(it, namefile, 'Delimiter', 'comma');

    % if plotFlag

    %     green = [0.4660 0.6740 0.1880];
    %     orange = [0.9290 0.6940 0.1250];

    %     figure;
    %     hold on
    %     set(0,'DefaultFigureWindowStyle','docked')

    %     p1 = semilogy(rrG(:,1), rrG(:,2), 'b-x','DisplayName',strcat('No deflation - PREC = ', num2str(PREC), ' - Restart = ', num2str(restart)),'linewidth', 2,'markersize', 10);
    %     p1.Color = green;
    %     if nbEigVec > 0
    %         p2 = semilogy(rrD(:,1), rrD(:,2), 'r-o','DisplayName',strcat('Deflation (', num2str(nbEigVec), ' vec) - PREC = ', num2str(PREC), ' - Restart = ', num2str(restart)),'linewidth', 2,'markersize', 10);
    %         p2.Color = orange;
    %     end
    %     set(gca, 'YScale', 'log')
    %     box on;
    %     grid on;
    %     ylim auto;
    %     xlabel('Iteration number', 'interpreter', 'latex', 'fontsize', 15);
    %     ylabel('Relative residual', 'interpreter', 'latex', 'fontsize', 15);
    %     title(['GMRES - "' benchmark '" - k = ' num2str(k) ' - h = ' num2str(h) ' - degree = ' num2str(degree) ' - PREC = ' num2str(PREC) ' - Dlambda = ' num2str(Dlambda) ' - Restart = ' num2str(restart)], 'interpreter', 'latex', 'fontsize', 20);
    %     legend('Location', 'southwest', 'fontsize', 15);
    % end

end