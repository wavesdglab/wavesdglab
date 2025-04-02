
clear;
clear global;

N=15;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);

global k h




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CAVITY BENCHMARK
k = 3.01*sqrt(2)*pi;
h = 1/32;
tol = 1e-6;
maxit = 2000;
L = 1;
degree = 2;
nbEigVec=1;
itout = 4;

run('cavity',degree,tol,maxit,itout,nbEigVec);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCATTERING OPEN CAVITY NEUMANN BENCHMARK
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
itout = 20;

run('scattering_openCavity_NEU',degree,tol,maxit,itout,nbEigVec);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%


function run(benchmark, degree, tol, maxit, itout, nbEigVec)

    global k h

    namediary = sprintf('output/log_mainHRV_%s_p%i_k=%g_def=%d.txt', benchmark, degree, k, nbEigVec);
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
    disp(['    Dlambda             ' num2str(Dlambda)]);
    disp(['---------------------------------------------------------']);

    [~, sysA] = computeSolNum2D_CG(mesh, dofm, 0);

    disp(['|  GMRES - No deflation']);
    [~, ~, itG, ~, ~, ~,hrv] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
    disp(['    converges in ' num2str(itG) ' iterations']);

    namefile = sprintf('output/hrv_%s_p%i_k=%g_def=0.csv', benchmark, degree, k);
    writematrix(hrv, namefile, 'Delimiter', 'comma');

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

        [Pdef,~,~] = computeDefOp(nbEigVec, eigvec, A);

        sysA.matA = Pdef*sysA.matA;
        sysA.rhsA = Pdef*sysA.rhsA;

        disp(['|  GMRES - Deflation = ' num2str(nbEigVec)]);
        [~, ~, itD, ~, ~, ~,hrvD] = solverGMRES_RP(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);
        disp(['    converges in ' num2str(itD) ' iterations']);

        namefile = sprintf('output/hrv_%s_p%i_k=%g_def=%d.csv', benchmark, degree, k, nbEigVec);
        writematrix(hrvD, namefile, 'Delimiter', 'comma');
    end


    diary off;
end

