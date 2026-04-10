%close all;
clear;
clear global;

N=15;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);

global Options

Options.Basis = 'Jacobi'; % Jacobi, Lagrange
Options.Error = 'L2'; % L2, H1

plotFlag = 1;
saveSolFlag = 1;
errorFlag = 1;

global omega cAir cObj rhoAir rhoObj h k Rdisk Rdom Rpml PML_TYPE

omega = 16.5962645;
cAir = 1.;
cObj = 2/3.;
rhoAir = 1.;
rhoObj = 1.;
h = 1/20;
tol = 1e-6;
degree = 2;

PML_TYPE = 'Circular';
Rdisk = 1.; Rdom = 1.2; Rpml = 4*h;

figure;
maxit = 5000; itout = 1;
PREC = 'none'; nbEigVec=[1]; restart = 0;
run('scattering_disk_penetrable',degree,PREC,tol,maxit,itout,nbEigVec,restart,plotFlag,saveSolFlag,errorFlag);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,PREC,tol,maxit,itout,nbEigVecList,restart,plotFlag,saveSolFlag,errorFlag)

global omega cAir cObj rhoAir rhoObj h k
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);

Dlambda = 2*pi/omega * (sqrt(dofm.numDofTRI) - 1);

disp(['---------------------------------------------------------']);
disp(['Method CG - Benchmark "' benchmark '"']);
disp(['---------------------------------------------------------']);
disp(['    omega               ' num2str(omega)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['    PREC                ' PREC]);
disp(['    Dlambda             ' num2str(Dlambda)]);
disp(['---------------------------------------------------------']);

switch PREC
    case {'none','eignum'}
        prec = 0;
    case {'CSL','ILU(CSL)','ILU(A)','geignum'}
        prec = 1;
    otherwise
        error('Error. \n%s is not a valid preconditioner', PREC);
end


[solA, sysA] = computeSolNum2D_CG_heterogeneous(mesh, dofm, PREC);


A = sysA.matA;
M = sysA.matP;
P = speye(size(A,2));
b = sysA.rhsA;

invPrec = @(x) P\x;
switch PREC
    case 'ILU(CSL)'
        [L,U] = ilu(P);
        invPrec = @(x) U\(L\x);
    case 'ILU(A)'
        [L,U] = ilu(A);
        invPrec = @(x) U\(L\x);
end

maxit = min(maxit, size(A,2));

if restart == 0
    m = size(A,2);
else
    m = restart;
    maxit = ceil(maxit/m);
end

for iterDeflation = 1:size(nbEigVecList(:),1)
    nbEigVec = nbEigVecList(iterDeflation);

    disp(['|  GMRES - PREC = ' PREC ' - Deflation = ' num2str(nbEigVec) ' - Restart = ' num2str(restart) '']);

    % Compute deflation matrices
    if nbEigVec > 0
        switch benchmark
            case 'cavity'
                [eigvec,nbEigVec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,'closestEigvec',omega);
            case 'scattering_openCavity_NEU'
                [eigvec,nbEigVec] = computeProjEigVec_openCavity_NEU(mesh, dofm, nbEigVec, omega);
            case 'scattering_openCavity_DIR'
                [eigvec,nbEigVec] = computeProjEigVec_openCavity_DIR(mesh, dofm, nbEigVec, omega);
            case 'scattering_disk_penetrable'
                [eigvec,nbEigVec] = computeProjEigVec_het(mesh, dofm, nbEigVec);
            otherwise
                error('Error. \n%s is not a valid benchmark for deflation', benchmark);
        end
%         switch PREC
%             case 'none'
%                 continue;
%             case 'eigvec'
%                 [eigvec,~] = eigs(A,nbEigVec,'sm');
%             case 'geignum'
%                 [eigvec,~] = eigs(A,M,nbEigVec,'sm');
%         end
    [PdefA,Pdef,Qdef,Q] = computeDefOp(nbEigVec, eigvec, A);
    [~,eigvalA] = eigs(A,5,'sm');
    [~,eigvalPdefA] = eigs(Pdef*A,5,'sm');
    else
        PdefA = @(x) A*x;
        Pdef = 1;
        Qdef = 1;
        Q = 0;
    end

    % Compute GMRES solution
    [uD, ~, ~, it, vecRes] = gmres(@(x) PdefA(invPrec(x)), Pdef*b, m, tol, maxit);
    disp(['    converged in ' num2str(size(vecRes,1)-1) ' iterations']);

    vecRes = vecRes(:)./vecRes(1);
    vecRes = vecRes(1:itout:end);
    vecIter = 0:itout:itout*size(vecRes,1)-1;
    xD = Q*b + Qdef*(invPrec(uD));

    if saveSolFlag
        namefile = sprintf('output/numSolD_%s_p%i_k=%g_prec=%s_def=%g_restart=%g.pos', benchmark, degree, omega, PREC, nbEigVec, restart);
        namesol = strcat('x_k=', num2str(omega), '_PREC=', PREC, '_def=', num2str(nbEigVec), '_restart=', num2str(restart));
        writeField2D(dofm, mesh, xD, namefile, namesol);
        system(['gmsh ../output/mesh.msh ' + namefile + '&']);
    end

    namefile = sprintf('output/historyGMRES_%s_p%i_k=%g_prec=%s_def=%g_restart=%g.csv', benchmark, degree, omega, PREC, nbEigVec, restart);
    writematrix([["it", "rrG"]; [vecIter' vecRes]], namefile, 'Delimiter', 'comma');

    if plotFlag
        hold on;
        semilogy(vecIter, vecRes, 'DisplayName',['PREC = ' PREC ' - nDef = ' num2str(nbEigVec) ' - nRestart = ' num2str(restart)]);
        set(gca, 'YScale', 'log')
        box on;
        grid on;
        ylim auto;
        xlabel('Iteration number');
        ylabel('Relative residual');
        title(['Benchmark "' benchmark '" - omega = ' num2str(omega) ' - h = ' num2str(h) ' - degree = ' num2str(degree) ' - Dlambda = ' num2str(Dlambda)]);
        legend('Location', 'southwest');
        ylim([0 maxit]);
        ylim([tol 1]);
    end

end

if errorFlag
    errorL2_iter = computeNormError2D_CG(mesh, dofm, xD);
    disp(['    L2-Error (iterative)   ' num2str(errorL2_iter,'%1.2e')]);
    errorL2_direct = computeNormError2D_CG(mesh, dofm, solA);
    disp(['    L2-Error (direct)   ' num2str(errorL2_direct,'%1.2e')]);
    solP = computeSolProjL2_2D_CG(mesh, dofm);
    errorProjL2 = computeNormError2D_CG(mesh, dofm, solP);
    disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
    disp(['---------------------------------------------------------']);

    writeField2D(dofm, mesh, xD, 'output/solIter.pos', "solIter");
    writeField2D(dofm, mesh, solA, 'output/solDirect.pos', "solDirect");
    writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
    writeField2D(dofm, mesh, xD-solA, 'output/errIter-Direct.pos', "errIterDirect");
    writeField2D(dofm, mesh, xD-solP, 'output/errIter-Ref.pos', "errIterRef");
    writeField2D(dofm, mesh, solA-solP, 'output/errDirect-Ref.pos', "errDirectRef");
    system('gmsh output/mesh.msh output/solRef.pos output/solDirect.pos output/errNum.pos output/errIter-Direct.pos output/errIter-Ref.pos output/errDirect-Ref.pos output/mesh.msh&');
end

% These variables must be cleared if 'cavity' is run after 'scattering'.
clear global LdomX LdomY LpmlX LpmlY;

end