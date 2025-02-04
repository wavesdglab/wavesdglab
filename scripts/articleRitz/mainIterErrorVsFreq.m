%close all;
clear;
clear global;

N=15;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);

global k h

plotFlag = 0;




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAVITY BENCHMARK
k = 3.01*sqrt(2)*pi;
h = 1/32;
rangeFreq = 12.7:0.009:13.8;
tol = 1e-10;
maxit = 2000;
L = 1;
degree = 2;
nbEigVec=1;


%% Without preconditioner
PREC = 'none';
restart = 0;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);


%% With preconditioner: CSL
PREC = 'CSL';
restart = 0;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 3;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 5;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);


%% With preconditioner: ILU(CSL)
PREC = 'ILU(CSL)';
restart = 0;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);


%% With preconditioner: ILU(A)
PREC = 'ILU(A)';
restart = 0;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('cavity',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCATTERING OPEN CAVITY DIRICHLET BENCHMARK
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
PREC = 'none';
restart = 0;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);


%% With preconditioner: CSL
PREC = 'CSL';
restart = 0;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 3;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 5;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);


%% With preconditioner: ILU(CSL)
PREC = 'ILU(CSL)';
restart = 0;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);


%% With preconditioner: ILU(A)
PREC = 'ILU(A)';
restart = 0;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('scattering_openCavity_DIR',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCATTERING OPEN CAVITY NEUMANN BENCHMARK
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
nbEigVec=1;


%% Without preconditioner
PREC = 'none';
restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);


%% With preconditioner: CSL
PREC = 'CSL';
restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 3;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 5;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);


%% With preconditioner: ILU(CSL)
PREC = 'ILU(CSL)';
restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);


%% With preconditioner: ILU(A)
PREC = 'ILU(A)';
restart = 0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 10;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);
restart = 25;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%


function run(benchmark,degree,PREC,tol,maxit,nbEigVec,restart,rangeFreq,plotFlag)

    global k h

    kmin = min(rangeFreq);
    kmax = max(rangeFreq);

    namediary = sprintf('output/log_iterErrorVsFreq_%s_p%i_range=%g-%g_prec=%s_def=%g_restart=%g.txt', benchmark, degree, kmin, kmax, PREC, nbEigVec, restart);
    delete(namediary);
    diary(namediary);

    if strcmp(benchmark, 'cavity')
        clear global LdomX LdomY LpmlX LpmlY
    end


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

    switch PREC
        case 'none'
            prec = 0;
        case 'CSL'
            prec = 1;
        case 'ILU(CSL)'
            prec = 1;
        case 'ILU(A)'
            prec = 1;
        otherwise
            error('Error. \n%s is not a valid preconditioner', PREC);
    end

    if nbEigVec > 0
        switch benchmark
            case 'cavity'
                [eigvec,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,'firstEigvec',k);

                errorL2 = zeros(length(rangeFreq),2);
                errorL2(:,1) = rangeFreq;

            case 'scattering_openCavity_NEU'
                [eigvec,nbEigVec] = computeProjEigVec_openCavity_NEU(mesh, dofm, nbEigVec, k);

            case 'scattering_openCavity_DIR'
                [eigvec,nbEigVec] = computeProjEigVec_openCavity_DIR(mesh, dofm, nbEigVec, k);

            otherwise
                error('Error. \n%s is not a valid benchmark for deflation', benchmark);

        end
    end


    for k = rangeFreq

        Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

        disp(['---------------------------------------------------------']);
        disp(['Method CG - Benchmark "' benchmark '"']);
        disp(['---------------------------------------------------------']);
        disp(['    k = ' num2str(k) ' - h = ' num2str(h) ' - Dlambda = ' num2str(Dlambda) ' - degree = ' num2str(degree) ' - PREC = ' PREC ' - Restart = ' num2str(restart)]);
        disp(['---------------------------------------------------------']);

        [~, sysA] = computeSolNum2D_CG(mesh, dofm, prec);

        A = sysA.matA;
        M = sysA.matP;
        b = sysA.rhsA;
        
        invPrec = @(x) M\x;

        if strcmp(PREC, 'ILU(CSL)')
            [L,U] = ilu(M);
            invPrec = @(x) U\(L\x);
        elseif strcmp(PREC, 'ILU(A)')
            [L,U] = ilu(A);
            invPrec = @(x) U\(L\x);
        end


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
        disp(['|  GMRES - No deflation - PREC = ' PREC ' - Restart = ' num2str(restart) '']);
        [uG, ~, ~, itG, ~] = gmres(@(x) A*(invPrec(x)),b,m,tol,maxit);
        itG = itG(2) + (itG(1)-1)*m;
        disp(['    converged in ' num2str(itG) ' iterations']);

        nbit(rangeFreq==k,2) = itG;

        if strcmp(benchmark, 'cavity')
            xG = invPrec(uG);
            errorL2(rangeFreq==k,2) = computeNormError2D_CG(mesh, dofm, xG);
        end

        if nbEigVec > 0
            % Deflation

            [Pdef,~,~] = computeDefOp(nbEigVec, eigvec, A);

            % Compute GMRES solution
            disp(['|  GMRES - Deflation - PREC = ' PREC ' - Restart = ' num2str(restart) '']);
            [~, ~, ~, itD, ~] = gmres(@(x) Pdef*(A*(invPrec(x))),Pdef*b,m,tol,maxit);
            itD = itD(2) + (itD(1)-1)*m;
            disp(['    converged in ' num2str(itD) ' iterations']);

            nbit(rangeFreq==k,3) = itD;
        end



    end

    namefile = sprintf('output/iterVsFreq_%s_p%i_range=%g-%g_prec=%s_def=%g_restart=%g.csv', benchmark, degree, kmin, kmax, PREC, nbEigVec, restart);
    writematrix([labels; nbit], namefile, 'Delimiter', 'comma');

    if strcmp(benchmark, 'cavity')
        namefile = sprintf('output/errorVsFreq_cavity_p%i_range=%g-%g_prec=%s_def=%g_restart=%g.csv', degree, kmin, kmax, PREC, nbEigVec, restart);
        writematrix([["k", "errorL2"]; errorL2], namefile, 'Delimiter', 'comma');
    end

    if plotFlag

        green = [0.4660 0.6740 0.1880];
        orange = [0.9290 0.6940 0.1250];

        figure;
        hold on;
        set(0,'DefaultFigureWindowStyle','docked');

        p1 = plot(nbit(:,1),nbit(:,2),'b-x','DisplayName',strcat('No deflation - PREC = ', PREC, ' - Restart = ', num2str(restart)),'linewidth', 2,'markersize', 10);
        p1.Color = green;
        if nbEigVec > 0
            p2 = plot(nbit(:,1),nbit(:,3),'r-o','DisplayName',strcat('Deflation (', num2str(nbEigVec), ' vec) - PREC = ', PREC, ' - Restart = ', num2str(restart)),'linewidth', 2,'markersize', 10);
            p2.Color = orange;
        end
        xlabel('k', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('Number of iterations', 'Interpreter', 'latex', 'FontSize', 15);
        title(['Nb of it vs freq - Benchmark "' benchmark '" - h = ' num2str(h) ' - degree = ' num2str(degree) ' - PREC = ' PREC ' - Restart = ' num2str(restart)], 'Interpreter', 'latex', 'FontSize', 20);
        legend('Location', 'best', 'FontSize', 15);
    end

    diary off;
end