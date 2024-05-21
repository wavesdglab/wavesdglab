
clear all;
%close all;

global k h

N=15;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);

computeSolNum2D = @computeSolNum2D_CG;

% Setup benchmark and parameters
benchmark = 'cavity';
switch benchmark
    case 'open'
        k = 15*pi;
        h = 1/16;
        tol = 1e-10; maxit = 1000; itout = 50;
    case 'cavity'
        k = 3.01*sqrt(2)*pi;
        h = 1/32;
        tol = 1e-6; maxit = 2000; itout =4;
        L = 1;
    case 'scatteringPML'
        k = 25;
        h = 0.05;
        tol = 1e-10; maxit = 2000; itout = 50;
        L = 1.1;
        R_disk = 1;
        L_PML = 0.2;
        computeSolNum2D = @computeSolNum2DPML_CG;
    case 'scattering_rec'
        global LdomX LdomY LpmlX LpmlY
        k = 10.5*pi;
        h = 1/8;
        tol = 1e-7; maxit = 5000; itout = 5;
        LdomX = 0.95;
        LdomY = 0.5;
        LpmlX = 0.8;
        LpmlY = 0.8;
    case 'waveguide'
        k = 6*pi;
        h = 1/8;
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 1; % P1
PREC = 0; % for preconditioner
% eigvecToDeflate = "closesteigvec"; %"firsteigvec" or "closesteigvec"
nbEigVec=1;


tabh = [1/8 1/16 1/32 1/64 1/128];
resvecGMRES = zeros(maxit, length(tabh));
% resvecADAnalyInterpol = zeros(maxit, length(tabh));
resvecADAnalyProj = zeros(maxit, length(tabh));
resvecADNum = zeros(maxit, length(tabh));
nbitGMRES = zeros(1, length(tabh));
% nbitADAnalyInterpol = zeros(1, length(tabh));
nbitADAnalyProj = zeros(1, length(tabh));
nbitADNum = zeros(1, length(tabh));

for h = tabh

    % Build mesh and DOF manager
    mesh = setupBenchmark2D(benchmark);
    mesh = buildConnectivity2D(mesh);
    dofm = buildDofManager2D_CG(mesh, degree); % espace fonctionnel discret

    Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1); % nb de points par longueur d'onde

    % -------------------------------------------------------------------------
    % Compute solution and error
    % -------------------------------------------------------------------------

    disp(['---------------------------------------------------------']);
    disp(['Method CG - Benchmark "' benchmark '"']);
    disp(['---------------------------------------------------------']);
    disp(['    k                   ' num2str(k)]);
    disp(['    h                   ' num2str(h)]);
    disp(['    degree              ' num2str(degree)]);
    disp(['    Dlambda             ' num2str(Dlambda)]);
    disp(['---------------------------------------------------------']);

    [~, sysA] = computeSolNum2D(mesh, dofm, PREC);

    A = sysA.matA;
    % M = sysA.matP;
    b = sysA.rhsA;

%     if maxit > size(A,2)
%         maxit = size(A,2);
%     end

    %%%%%%%%%%% No deflation %%%%%%%%%%%

    % Compute GMRES without prec
    [xGMRES, ~, ~, itGMRES, rrGMRES] = gmres(A,b,[],tol,maxit);
    itGMRES = itGMRES(2);
    rrGMRES = rrGMRES(:)./rrGMRES(1);
    rrGMRES = rrGMRES(1:itout:end);

    resvecGMRES(1:length(rrGMRES), tabh == h) = rrGMRES;
    nbitGMRES(tabh == h) = itGMRES;



    %%%%%%%%%%% Deflation analytique %%%%%%%%%%%

    [eigvec,nbEigVec] = computeEigVec2D_cavity(mesh, dofm, nbEigVec,'closesteigvec');
    % [eigvec,~] = eigs(A,nbEigVec,'smallestabs');

    [P,Q] = computeDefOp(nbEigVec, eigvec, A);

    % [~, evdef] = eigs((P+Q)*A,50-nbEigVec,'smallestabs');
    % evdef = diag(evdef);



    %%% No preconditioner :

    % Compute GMRES with DEF1 and closest eigvec : P*A*x = P*b
    % [xD, ~, ~, itD, rrD] = gmres(P*A,P*b,[],tol,maxit);
    % itD = itD(2);
    % rrD = rrD(:)./rrD(1);
    % rrD = rrD(1:itout:end);
    % xD = P'*xD + Q*b;



    % Compute GMRES with ADEF1 and closest eigvec : (P+Q)*A*x = (P+Q)*b
    [xAD, ~, ~, itAD, rrAD] = gmres((P+Q)*A,(P+Q)*b,[],tol,maxit);
    itAD = itAD(2);
    rrAD = rrAD(:)./rrAD(1);
    rrAD = rrAD(1:itout:end);

    resvecADAnalyProj(1:length(rrAD), tabh == h) = rrAD;
    nbitADAnalyProj(tabh == h) = itAD;

    writeField2D(dofm, mesh, eigvec(:,1), "output/eigenvec_analy" + string(h) + ".pos", "eigenvec_analy" + string(h));

    %%%%%%%%%%% Deflation numérique %%%%%%%%%%%

    [eigvec,~] = eigs(A,nbEigVec,'smallestabs');
    [P,Q] = computeDefOp(nbEigVec, eigvec, A);
    % Compute GMRES with ADEF1 and closest eigvec : (P+Q)*A*x = (P+Q)*b
    [xAD, ~, ~, itAD, rrAD] = gmres((P+Q)*A,(P+Q)*b,[],tol,maxit);
    itAD = itAD(2);
    rrAD = rrAD(:)./rrAD(1);
    rrAD = rrAD(1:itout:end);

    resvecADNum(1:length(rrAD), tabh == h) = rrAD;
    nbitADNum(tabh == h) = itAD;

    disp(['---------------------------------------------------------']);
    disp(['Number of iterations GMRES: ' num2str(itGMRES)]);
    disp(['Number of iterations GMRES with ADEF analytique: ' num2str(nbitADAnalyProj(tabh == h))]);
    disp(['Number of iterations GMRES with ADEF numérique: ' num2str(nbitADNum(tabh == h))]);
    disp(['---------------------------------------------------------']);

    writeField2D(dofm, mesh, eigvec(:,1), "output/eigenvec_num" + string(h) + ".pos", "eigenvec_num" + string(h));

    maxIt = max(max(nbitGMRES), max(max(nbitADAnalyProj), max(nbitADNum)));
    minIt = min(min(nbitGMRES), min(min(nbitADAnalyProj), min(nbitADNum)));

    figure
    hold on
    set(0,'DefaultFigureWindowStyle','docked')

    p1 = semilogy(1:nbitGMRES(tabh == h),resvecGMRES(1:nbitGMRES(tabh == h),tabh == h),'b-o','DisplayName','rrGMRES','linewidth', 2,'markersize', 10);
    p2 = semilogy(1:nbitADAnalyProj(tabh == h),resvecADAnalyProj(1:nbitADAnalyProj(tabh == h),tabh == h),'r-+','DisplayName','rrADAnalyProj','linewidth', 2,'markersize', 10);
    p3 = semilogy(1:nbitADNum(tabh == h),resvecADNum(1:nbitADNum(tabh == h),tabh == h),'g-x','DisplayName','rrADNum','linewidth', 2,'markersize', 10);

    set(gca, 'YScale', 'log')
    box on
    grid on
    xlim([0 maxIt+1]);
    ylim auto;
    title(['CG - ' benchmark ' - GMRES - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree) ' - nbEigvec=' num2str(nbEigVec)], 'interpreter', 'latex', 'fontsize', 20)
    xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
    ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
    legend('Location', 'southwest', 'fontsize', 15)

end

%%% Plot results %%%

% green = [0.4660 0.6740 0.1880];
% magenta = [0.4940 0.1840 0.5560];
% orange = [0.9290 0.6940 0.1250];
% cyan = [0.3010 0.7450 0.9330];




%%% Residuals

% maxIt = max(max(nbitGMRES), max(max(nbitADAnalyProj), max(nbitADNum)));
% minIt = min(min(nbitGMRES), min(min(nbitADAnalyProj), min(nbitADNum)));


% figure
% hold on
% set(0,'DefaultFigureWindowStyle','docked')

% p1 = semilogy(1:nbitGMRES(1),resvecGMRES(1:nbitGMRES(1),1),'b-o','DisplayName','rrGMRES h=1/8','linewidth', 2,'markersize', 10);
% p2 = semilogy(1:nbitADAnalyProj(1),resvecADAnalyProj(1:nbitADAnalyProj(1),1),'b-+','DisplayName','rrADAnalyProj h=1/8','linewidth', 2,'markersize', 10);
% p3 = semilogy(1:nbitADNum(1),resvecADNum(1:nbitADNum(1),1),'b-x','DisplayName','rrADNum h=1/8','linewidth', 2,'markersize', 10);

% p4 = semilogy(1:nbitGMRES(2),resvecGMRES(1:nbitGMRES(2),2),'r-o','DisplayName','rrGMRES h=1/16','linewidth', 2,'markersize', 10);
% p5 = semilogy(1:nbitADAnalyProj(2),resvecADAnalyProj(1:nbitADAnalyProj(2),2),'r-+','DisplayName','rrADAnalyProj h=1/16','linewidth', 2,'markersize', 10);
% p6 = semilogy(1:nbitADNum(2),resvecADNum(1:nbitADNum(2),2),'r-x','DisplayName','rrADNum h=1/16','linewidth', 2,'markersize', 10);

% p7 = semilogy(1:nbitGMRES(3),resvecGMRES(1:nbitGMRES(3),3),'g-o','DisplayName','rrGMRES h=1/32','linewidth', 2,'markersize', 10);
% % p7.Color = green;
% p8 = semilogy(1:nbitADAnalyProj(3),resvecADAnalyProj(1:nbitADAnalyProj(3),3),'g-+','DisplayName','rrADAnalyProj h=1/32','linewidth', 2,'markersize', 10);
% % p8.Color = green;
% p9 = semilogy(1:nbitADNum(3),resvecADNum(1:nbitADNum(3),3),'g-x','DisplayName','rrADNum h=1/32','linewidth', 2,'markersize', 10);
% % p9.Color = green;

% p10 = semilogy(1:nbitGMRES(4),resvecGMRES(1:nbitGMRES(4),4),'m-o','DisplayName','rrGMRES h=1/64','linewidth', 2,'markersize', 10);
% % p10.Color = magenta;
% p11 = semilogy(1:nbitADAnalyProj(4),resvecADAnalyProj(1:nbitADAnalyProj(4),4),'m-+','DisplayName','rrADAnalyProj h=1/64','linewidth', 2,'markersize', 10);
% % p11.Color = magenta;
% p12 = semilogy(1:nbitADNum(4),resvecADNum(1:nbitADNum(4),4),'m-x','DisplayName','rrADNum h=1/64','linewidth', 2,'markersize', 10);
% % p12.Color = magenta;

% p13 = semilogy(1:nbitGMRES(5),resvecGMRES(1:nbitGMRES(5),5),'c-o','DisplayName','rrGMRES h=1/128','linewidth', 2,'markersize', 10);
% % p13.Color = cyan;
% p14 = semilogy(1:nbitADAnalyProj(5),resvecADAnalyProj(1:nbitADAnalyProj(5),5),'c-+','DisplayName','rrADAnalyProj h=1/128','linewidth', 2,'markersize', 10);
% % p14.Color = cyan;
% p15 = semilogy(1:nbitADNum(5),resvecADNum(1:nbitADNum(5),5),'c-x','DisplayName','rrADNum h=1/128','linewidth', 2,'markersize', 10);
% % p15.Color = cyan;

% set(gca, 'YScale', 'log')
% box on
% grid on
% xlim([0 maxIt+1]);
% ylim auto;
% title(['CG - ' benchmark ' - GMRES - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree) ' - nbEigvec=' num2str(nbEigVec)], 'interpreter', 'latex', 'fontsize', 20)
% xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
% ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
% legend('Location', 'southwest', 'fontsize', 15)