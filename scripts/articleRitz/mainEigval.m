%close all;
clear;
clear global;

N=32;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);

global k h

plotFlag = 1;
saveSolFlag = 1;




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAVITY BENCHMARK
k = 3.01*sqrt(2)*pi;
h = 1/32;
L = 1;
degree = 2;
nbEigVec=1;


%% Without preconditioner
PREC = 'none';
run('cavity',degree,PREC,nbEigVec,plotFlag);

% %% With preconditioner: CSL
% PREC = 'CSL';
% run('cavity',degree,PREC,nbEigVec,plotFlag);


%% With preconditioner: ILU(CSL)
PREC = 'ILU(CSL)';
run('cavity',degree,PREC,nbEigVec,plotFlag);

% %% With preconditioner: ILU(A)
% PREC = 'ILU(A)';
% run('cavity',degree,PREC,nbEigVec,plotFlag);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCATTERING OPEN CAVITY NEUMANN BENCHMARK
global LdomX LdomY LpmlX LpmlY
k = 23.591;
% k = 23.82;
% k = 24.275;
h = 1/20;
LdomX = 0.95;
LdomY = 0.5;
LpmlX = 0.2;
LpmlY = 0.2;
degree = 3;
nbEigVec=1;


%% Without preconditioner
PREC = 'none';
run('scattering_openCavity_NEU',degree,PREC,nbEigVec,plotFlag);


% %% With preconditioner: CSL
% PREC = 'CSL';
% run('scattering_openCavity_NEU',degree,PREC,nbEigVec,plotFlag);


%% With preconditioner: ILU(CSL)
PREC = 'ILU(CSL)';
run('scattering_openCavity_NEU',degree,PREC,nbEigVec,plotFlag);

% %% With preconditioner: ILU(A)
% PREC = 'ILU(A)';
% run('scattering_openCavity_NEU',degree,PREC,nbEigVec,plotFlag);




% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCATTERING OPEN CAVITY DIRICHLET BENCHMARK
% global LdomX LdomY LpmlX LpmlY WRITE_FIELD_ABSOLUTE
% WRITE_FIELD_ABSOLUTE = 1;
% k = 23.676;
% h = 1/20;
% LdomX = 0.95;
% LdomY = 0.5;
% LpmlX = 0.2;
% LpmlY = 0.2;
% degree = 3;
% nbEigVec = 1;


% %% Without preconditioner
% PREC = 'none';
% run('scattering_openCavity_DIR',degree,PREC,nbEigVec,plotFlag);


% %% With preconditioner: CSL
% PREC = 'CSL';
% run('scattering_openCavity_DIR',degree,PREC,nbEigVec,plotFlag);


% %% With preconditioner: ILU(CSL)
% PREC = 'ILU(CSL)';
% run('scattering_openCavity_DIR',degree,PREC,nbEigVec,plotFlag);


% %% With preconditioner: ILU(A)
% PREC = 'ILU(A)';
% run('scattering_openCavity_DIR',degree,PREC,nbEigVec,plotFlag);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%


function run(benchmark,degree,PREC,nbEigVec,plotFlag)

    global k h

    namediary = sprintf('output/log_mainEigval_%s_p%i_k=%g_prec=%s_def=%g.txt', benchmark, degree, k, PREC, nbEigVec);
    if exist(namediary, 'file')
        delete(namediary);
    end
    diary(namediary);

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
        case 'ILU(CSL)'
            prec = 1;
        case 'ILU(A)'
            prec = 1;
        otherwise
            error('Error. \n%s is not a valid preconditioner', PREC);
    end

    [~, sysA] = computeSolNum2D_CG(mesh, dofm, prec);

    A = sysA.matA;
    M = sysA.matP;

    [~,eigvalA] = eigs(A/M,size(A,2),'sm');
    eigvalA = diag(eigvalA);
    eigvalA = sort(eigvalA,'descend');
    eigvalA = eigvalA(1:end);

    writematrix(eigvalA, sprintf('output/eigvalA_%s_p%i_k=%g_prec=%s_def=%g.csv', benchmark, degree, k, PREC, nbEigVec), 'Delimiter', 'comma');


    if plotFlag
        figure;
        plot(real(eigvalA),imag(eigvalA),'o');
        xlabel('Re');
        ylabel('Im');
        title('Eigenvalues of the matrix A');
    end


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

        [~,Pdef,~,~] = computeDefOp(nbEigVec, eigvec, A);


        % Deflation

        [~,eigvalPA] = eigs(Pdef*A/M,size(A,2),'sm');
        eigvalPA = diag(eigvalPA);
        eigvalPA = sort(eigvalPA,'descend');
        eigvalPA = eigvalPA(1:end);

        writematrix(eigvalPA, sprintf('output/eigvalPA_%s_p%i_k=%g_prec=%s_def=%g.csv', benchmark, degree, k, PREC, nbEigVec), 'Delimiter', 'comma');

        if plotFlag
            figure;
            plot(real(eigvalPA),imag(eigvalPA),'o');
            xlabel('Re');
            ylabel('Im');
            title('Eigenvalues of the matrix P*A');
        end

    end

    diary off;
end
