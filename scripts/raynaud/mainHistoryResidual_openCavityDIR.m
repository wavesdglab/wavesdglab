
clear all;
%close all;

global k h

N=15;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);


% Setup benchmark and parameters
benchmark = 'scattering_openCavity_DIR';
global LdomX LdomY LpmlX LpmlY
k = 23.676;
tol = 1e-6;
maxit = 5000;
itout = 10;
LdomX = 0.95;
LdomY = 0.5;
LpmlX = 0.2;
LpmlY = 0.2;
degree = 3;
PREC = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = 1/8;
nbEigVec=1;
run(benchmark,degree,PREC,tol,maxit,itout,nbEigVec);

h = 1/16;
nbEigVec=1;
run(benchmark,degree,PREC,tol,maxit,itout,nbEigVec);

h = 1/32;
nbEigVec=1;
run(benchmark,degree,PREC,tol,maxit,itout,nbEigVec);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function run(benchmark,degree,PREC,tol,maxit,itout,nbEigVec)
    global k h

    % Build mesh and DOF manager
    mesh = setupBenchmark2D(benchmark);
    mesh = buildConnectivity2D(mesh);
    dofm = buildDofManager2D_CG(mesh, degree);

    Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

    % -------------------------------------------------------------------------
    % Compute solution
    % -------------------------------------------------------------------------

    disp(['---------------------------------------------------------']);
    disp(['Method CG - Benchmark "' benchmark '"']);
    disp(['---------------------------------------------------------']);
    disp(['    k                   ' num2str(k)]);
    disp(['    h                   ' num2str(h)]);
    disp(['    degree              ' num2str(degree)]);
    disp(['    Dlambda             ' num2str(Dlambda)]);
    disp(['---------------------------------------------------------']);

    disp(['| Compute system and deflation subspace with analytical eigenvectors...']);

    [~, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);

    A = sysA.matA;
    M = sysA.matP;
    b = sysA.rhsA;
    AMinv = A/M;

    [eigvec,nbEigVec] = computeProjEigVec_openCavity(mesh, dofm, nbEigVec, k);
    [P,Q] = computeDefOp(nbEigVec, eigvec, A);

    MinvP = M\P;

    if maxit > size(A,2)
        maxit = size(A,2);
    end


    disp(['|             Done']);

    % disp(['| Compute eigenvalues of A...']);
    % [~, evA] = eigs(A,size(A,2),'smallestabs');
    % evA = diag(evA);
    % evA = [real(evA) imag(evA)];
    % csvwrite(["output/evA.csv"],evA);

    % disp(['|             Done']);


    % disp(['| Compute eigenvalues of A/M...']);
    % [~, evAMinv] = eigs(AMinv,50,'smallestabs');
    % evAMinv = diag(evAMinv);
    % evAMinv = [real(evAMinv) imag(evAMinv)];
    % csvwrite(["output/spectrum_prec.csv"],evAMinv);
    % disp(['|             Done']);


    % disp(['| Compute eigenvalues of A*(P+Q)...']);
    % [~, evD] = eigs(A*(P+Q),50-nbEigVec,'smallestabs');
    % evD = diag(evD);
    % evD = [real(evD) imag(evD)];
    % disp(['|             Done']);

    % disp(['| Compute eigenvalues of A*(M\P+Q)...']);
    % [~, evPD] = eigs(A*(MinvP+Q),50-nbEigVec,'smallestabs');
    % evPD = diag(evPD);
    % evPD = [real(evPD) imag(evPD)];
    % disp(['|             Done']);


    %%%%%%%%%%% No deflation %%%%%%%%%%%

    % Compute GMRES with prec
    disp(['| Preconditioned GMRES...']);
    [xGMRESP, ~, ~, itGMRESP, rrGMRESP] = gmres(AMinv,b,[],tol,maxit);
    itGMRESP = itGMRESP(2);
    rrGMRESP = rrGMRESP(:)./rrGMRESP(1);
    rrGMRESP = rrGMRESP(1:itout:end);
    iterGMRESP = 0:itout:itout*size(rrGMRESP,1)-1;
    rrGMRESP = [iterGMRESP' rrGMRESP];
    xGMRESP = M\xGMRESP;
    disp(['|             converges in ' num2str(itGMRESP) ' iterations']);



    % Compute GMRES without prec
    disp(['| GMRES...']);
    [xGMRES, ~, ~, itGMRES, rrGMRES] = gmres(A,b,[],tol,maxit);
    itGMRES = itGMRES(2);
    rrGMRES = rrGMRES(:)./rrGMRES(1);
    rrGMRES = rrGMRES(1:itout:end);
    iterGMRES = 0:itout:itout*size(rrGMRES,1)-1;
    rrGMRES = [iterGMRES' rrGMRES];
    disp(['|             converges in ' num2str(itGMRES) ' iterations']);



    %%%%%%%%%%% Deflation %%%%%%%%%%%


    %%% No preconditioner :


    % Compute GMRES with ADEF1 and closest eigvec : A*(P+Q)*u = b, x = (P+Q)*u
    disp(['| GMRES with ADEF1...']);
    [uADana, ~, ~, itADana, rrADana] = gmres(A*(P+Q),b,[],tol,maxit);
    itADana = itADana(2);
    rrADana = rrADana(:)./rrADana(1);
    rrADana = rrADana(1:itout:end);
    iterAD = 0:itout:itout*size(rrADana,1)-1;
    rrADana = [iterAD' rrADana];
    xADana = (P+Q)*uADana;
    disp(['|             converges in ' num2str(itADana) ' iterations']);




    %%% Add preconditioner :



    % Compute GMRES with ADEF1 and closest eigvec and prec : A*(M\P+Q)*u = b, x = (M\P+Q)*u
    disp(['| Preconditoned  GMRES with ADEF1...']);
    [uPADana, ~, ~, itPADana, rrPADana] = gmres(A*(MinvP+Q),b,[],tol,maxit);
    itPADana = itPADana(2);
    rrPADana = rrPADana(:)./rrPADana(1);
    rrPADana = rrPADana(1:itout:end);
    iterPAD = 0:itout:itout*size(rrPADana,1)-1;
    rrPADana = [iterPAD' rrPADana];
    xPADana = (MinvP+Q)*uPADana;
    disp(['|             converges in ' num2str(itPADana) ' iterations']);


    disp(['| Compute deflation subspace with numerical eigenvectors...']);

    [eigvec,~] = eigs(A,nbEigVec,'smallestabs');
    [P,Q] = computeDefOp(nbEigVec, eigvec, A);

    MinvP = M\P;

    disp(['|             Done']);

    %%%%%%%%%%% Deflation %%%%%%%%%%%

    % Compute GMRES with ADEF1 and closest eigvec : A*(P+Q)*u = b, x = (P+Q)*u
    disp(['| GMRES with ADEF1 and numerical eigenvectors...']);
    [uADnum, ~, ~, itADnum, rrADnum] = gmres(A*(P+Q),b,[],tol,maxit);
    itADnum = itADnum(2);
    rrADnum = rrADnum(:)./rrADnum(1);
    rrADnum = rrADnum(1:itout:end);
    iterADnum = 0:itout:itout*size(rrADnum,1)-1;
    rrADnum = [iterADnum' rrADnum];
    xADnum = (P+Q)*uADnum;
    disp(['|             converges in ' num2str(itADnum) ' iterations']);

    % Compute GMRES with ADEF1 and closest eigvec and prec : A*(M\P+Q)*u = b, x = (M\P+Q)*u
    disp(['| Preconditioned GMRES with ADEF1 and numerical eigenvectors...']);
    [uPADnum, ~, ~, itPADnum, rrPADnum] = gmres(A*(MinvP+Q),b,[],tol,maxit);
    itPADnum = itPADnum(2);
    rrPADnum = rrPADnum(:)./rrPADnum(1);
    rrPADnum = rrPADnum(1:itout:end);
    iterPADnum = 0:itout:itout*size(rrPADnum,1)-1;
    rrPADnum = [iterPADnum' rrPADnum];
    xPADnum = (MinvP+Q)*uPADnum;
    disp(['|             converges in ' num2str(itPADnum) ' iterations']);

    % compute incident field
    disp(['| Compute incident field...']);

    solInc = -computeSolProjL2_2D_CG(mesh, dofm);

    %%% Save results %%%

    disp(['| Save results...']);

    it = [itGMRES itGMRESP itADana itPADana itADnum itPADnum];

    folder = "output/mainHistoryResidual_openCavityDIR";
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    global WRITE_FIELD_ABSOLUTE
    WRITE_FIELD_ABSOLUTE = 1;

    % for i=1:nbEigVec
    %     writeField2D(dofm, mesh, eigvec(:,i), folder+"/eigvec"+num2str(i)+".pos", "eigvec"+num2str(i));
    % end

    writeField2D(dofm, mesh, xGMRESP+solInc, folder+"/solGMRESP.pos", "solGMRESP");
    writeField2D(dofm, mesh, xGMRES+solInc, folder+"/solGMRES.pos", "solGMRES");
    writeField2D(dofm, mesh, xADana+solInc, folder+"/solADefana.pos", "solADefana");
    writeField2D(dofm, mesh, xPADana+solInc, folder+"/solADefPana.pos", "solADefPana");
    writeField2D(dofm, mesh, xADnum+solInc, folder+"/solADefnum.pos", "solADefnum");
    writeField2D(dofm, mesh, xPADnum+solInc, folder+"/solADefPnum.pos", "solADefPnum");

    csvwrite([folder+"/rrGMRES.csv"],rrGMRES);
    csvwrite([folder+"/rrGMRESP.csv"],rrGMRESP);
    csvwrite([folder+"/rrADana.csv"],rrADana);
    csvwrite([folder+"/rrPADana.csv"],rrPADana);
    csvwrite([folder+"/rrADnum.csv"],rrADnum);
    csvwrite([folder+"/rrPADnum.csv"],rrPADnum);


    % csvwrite([folder+"/evA.csv"],evA);
    % csvwrite([folder+"/evAMinv.csv"],evAMinv);
    % csvwrite([folder+"/evD.csv"],evD);
    % csvwrite([folder+"/evPD.csv"],evPD);

    % csvwrite([folder+"/it.csv"],it');

    %%% Plot results %%%

    green = [0.4660 0.6740 0.1880];
    magenta = [0.4940 0.1840 0.5560];
    orange = [0.9290 0.6940 0.1250];
    cyan = [0.3010 0.7450 0.9330];

    %%% Spectra


    % figure;
    % hold on;
    % s1 = scatter(real(evA),imag(evA),100,'DisplayName','Eigenvalues of A');
    % s1.Marker = '+';
    % s1.MarkerEdgeColor = 'b';
    %
    % s2 = scatter(real(evAMinv),imag(evAMinv),100, 'DisplayName','Eigenvalues of A/M');
    % s2.Marker = 'x';
    % s2.MarkerEdgeColor = 'k';
    %
    % s3 = scatter(real(evD),imag(evD),100, 'DisplayName','Eigenvalues of A*(P+Q)');
    % s3.Marker = 'o';
    % s3.MarkerEdgeColor = 'r';
    %
    % s4 = scatter(real(evPD),imag(evPD),100, 'DisplayName','Eigenvalues of A*(M\P+Q)');
    % s4.Marker = 'o';
    % s4.MarkerEdgeColor = 'g';
    %
    % legend('Location', 'southwest', 'fontsize', 15)
    %
    % grid on; box on;
    % title(['Spectra of $A$ - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree)], 'interpreter', 'latex', 'fontsize', 20)
    % % xlim([-0.04 0.06]);
    % % hold on; plot(fovals(A));
    % set(0,'DefaultFigureWindowStyle','docked')

    disp(['|             Done']);

    %%% Residuals

    maxIt = max(it);
    minIt = min(it);


    figure
    hold on
    set(0,'DefaultFigureWindowStyle','docked')

    p1 = semilogy(rrGMRES(:,1),rrGMRES(:,2),'b-o','DisplayName','Relative residual','linewidth', 2,'markersize', 10);
    p2 = semilogy(rrGMRESP(:,1),rrGMRESP(:,2),'r-o','DisplayName','Relative residual with CSL','linewidth', 2,'markersize', 10);
    p3 = semilogy(rrADana(:,1),rrADana(:,2),'c-o','DisplayName','Relative residual with ADEF1','linewidth', 2,'markersize', 10);
    p3.Color = cyan;
    p4  = semilogy(rrPADana(:,1),rrPADana(:,2),'y-o','DisplayName','Relative residual with ADEF1 and CSL','linewidth', 2,'markersize', 10);
    p4.Color = orange;
    p5 = semilogy(rrADnum(:,1),rrADnum(:,2),'c-x','DisplayName','Relative residual with ADEF1 and numerical eigenvectors','linewidth', 2,'markersize', 10);
    p5.Color = cyan;
    p6 = semilogy(rrPADnum(:,1),rrPADnum(:,2),'y-x','DisplayName','Relative residual with ADEF1 and CSL and numerical eigenvectors','linewidth', 2,'markersize', 10);
    p6.Color = orange;

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