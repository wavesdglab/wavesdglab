%close all;
clear;
clear global;

global k h

plotFlag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Benchmark 'Cavity'
% Alternatives for PREC: 'none', 'CSL', 'ILU(CSL)', 'ILU(A)'

h = 1/32;
rangeFreq = 12.7:0.01:13.8;
tol = 1e-6;
maxit = 2000;
degree = 2;
restart = 0;

figure;
PREC = 'none'; defFreq = 3.01*sqrt(2)*pi; nbEigVec = 0;
run('cavity',degree,PREC,tol,maxit,nbEigVec,defFreq,restart,rangeFreq,plotFlag);
PREC = 'none'; defFreq = 3.01*sqrt(2)*pi; nbEigVec = 1;
run('cavity',degree,PREC,tol,maxit,nbEigVec,defFreq,restart,rangeFreq,plotFlag);

figure;
PREC = 'ILU(CSL)'; defFreq = 3.01*sqrt(2)*pi; nbEigVec = 0;
run('cavity',degree,PREC,tol,maxit,nbEigVec,defFreq,restart,rangeFreq,plotFlag);
PREC = 'ILU(CSL)'; defFreq = 3.01*sqrt(2)*pi; nbEigVec = 1;
run('cavity',degree,PREC,tol,maxit,nbEigVec,defFreq,restart,rangeFreq,plotFlag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benchmark 'Scattering Open Cavity Neumann'
% Alternative with Dirichlet: 'scattering_openCavity_DIR' with 23.676

global LdomX LdomY LpmlX LpmlY
LdomX = 0.95; LdomY = 0.5; LpmlX = 0.2; LpmlY = 0.2;

h = 1/20;
rangeFreq = 23.5:0.1:24.5;
tol = 1e-6;
maxit = 2000;
degree = 3;
restart = 0;

figure;
PREC = 'none'; defFreq = 23.591; nbEigVec=0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,defFreq,restart,rangeFreq,plotFlag);
PREC = 'none'; defFreq = 23.591; nbEigVec=1;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,defFreq,restart,rangeFreq,plotFlag);
return;
figure;
PREC = 'ILU(CSL)'; defFreq = 23.591; nbEigVec=0;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,defFreq,restart,rangeFreq,plotFlag);
PREC = 'ILU(CSL)'; defFreq = 23.591; nbEigVec=1;
run('scattering_openCavity_NEU',degree,PREC,tol,maxit,nbEigVec,defFreq,restart,rangeFreq,plotFlag);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,PREC,tol,maxit,nbEigVec,defFreq,restart,rangeFreq,plotFlag)

global k h

kmin = min(rangeFreq);
kmax = max(rangeFreq);

mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

nbit = zeros(length(rangeFreq),2);
labels = ["k", "it"];
nbit(:,1) = rangeFreq;

if strcmp(benchmark, 'cavity')
    errorL2 = zeros(length(rangeFreq),2);
    errorL2(:,1) = rangeFreq;
end

switch PREC
    case 'none'
        prec = 0;
    case {'CSL', 'ILU(CSL)', 'ILU(A)'}
        prec = 1;
    otherwise
        error('Error. \n%s is not a valid preconditioner', PREC);
end

if nbEigVec > 0
    switch benchmark
        case 'cavity'
            [eigvec,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,'closestEigvec',defFreq);
        case 'scattering_openCavity_NEU'
            [eigvec,nbEigVec] = computeProjEigVec_openCavity_NEU(mesh, dofm, nbEigVec, defFreq);
        case 'scattering_openCavity_DIR'
            [eigvec,nbEigVec] = computeProjEigVec_openCavity_DIR(mesh, dofm, nbEigVec, defFreq);
        otherwise
            error('Error. \n%s is not a valid benchmark for deflation', benchmark);
            
    end
end

disp(['---------------------------------------------------------']);
disp(['Method CG - Benchmark "' benchmark '"']);

for k = rangeFreq
    
    Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);
    
    disp(['---------------------------------------------------------']);
    disp(['    k = ' num2str(k) ' - h = ' num2str(h) ' - Dlambda = ' num2str(Dlambda) ' - degree = ' num2str(degree) ' - PREC = ' PREC ' - Restart = ' num2str(restart)]);
    
    [~, sysA] = computeSolNum2D_CG(mesh, dofm, prec);
    
    A = sysA.matA;
    M = sysA.matP;
    b = sysA.rhsA;
    
    invPrec = @(x) M\x;
    
    switch PREC
        case 'ILU(CSL)'
            [L,U] = ilu(M);
            invPrec = @(x) U\(L\x);
        case 'ILU(A)'
            [L,U] = ilu(A);
            invPrec = @(x) U\(L\x);
    end
    
    
    if k == rangeFreq(1)
        maxit = min(maxit, size(A,2));
        if restart == 0
            m = size(A,2);
        else
            m = restart;
            maxit = ceil(maxit/m);
        end
    end
    
    if nbEigVec > 0
        [PdefA,Pdef,~,~] = computeDefOp(nbEigVec, eigvec, A);
    else
        PdefA = @(x) A*x;
        Pdef = speye(size(A,1));
    end
    
    % Compute GMRES solution
    [u, ~, ~, it, ~] = gmres(@(x) PdefA(invPrec(x)),Pdef*b,m,tol,maxit);
    it = it(2) + (it(1)-1)*m;
    disp(['|  GMRES converged in ' num2str(it) ' iterations']);
    
    nbit(rangeFreq==k,2) = it;
    
    if strcmp(benchmark, 'cavity')
        x = invPrec(u);
        errorL2(rangeFreq==k,2) = computeNormError2D_CG(mesh, dofm, x);
    end
    
    
end

namefile = sprintf('output/iterVsFreq_%s_p%i_range=%g-%g_prec=%s_def=%g_restart=%g.csv', benchmark, degree, kmin, kmax, PREC, nbEigVec, restart);
writematrix([labels; nbit], namefile, 'Delimiter', 'comma');

if strcmp(benchmark, 'cavity')
    namefile = sprintf('output/errorVsFreq_cavity_p%i_range=%g-%g_prec=%s_def=%g_restart=%g.csv', degree, kmin, kmax, PREC, nbEigVec, restart);
    writematrix([["k", "errorL2"]; errorL2], namefile, 'Delimiter', 'comma');
end

if plotFlag
    hold on;
    
    plot(nbit(:,1),nbit(:,2),'b-x','DisplayName',strcat('ndef = ', num2str(nbEigVec), ' - PREC = ', PREC, ' - Restart = ', num2str(restart)),'linewidth', 2,'markersize', 10);
    box on;
    grid on;
    ylim auto;
    xlabel('k', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Number of iterations', 'Interpreter', 'latex', 'FontSize', 15);
    title(['Nb of it vs freq - Benchmark "' benchmark '" - h = ' num2str(h) ' - degree = ' num2str(degree) ' - PREC = ' PREC ' - Restart = ' num2str(restart)], 'Interpreter', 'latex', 'FontSize', 20);
    legend('Location', 'best', 'FontSize', 15);
end

% These variables must be cleared if 'cavity' is run after 'scattering'.
clear global LdomX LdomY LpmlX LpmlY;
end