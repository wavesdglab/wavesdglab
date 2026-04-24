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

plotFlag = 0;
generalizedFlag = 1;
generalizedPMLFlag = 1;

global omega cAir cObj rhoAir rhoObj h k Rdisk Rdom Rpml PML_TYPE

omega = 16.5962645;
% omega = 6.2;
cAir = 1.;
cObj = 2/3.;
rhoAir = 1.;
rhoObj = 1.;
% h = 1/25;
h = 1/10;
tol = 1e-6;
degree = 2;

PML_TYPE = 'Circular';
Rdisk = 1.; Rdom = 1.2; Rpml = 4*h;

figure;
maxit = 5000;
PREC = 'none'; nbEigVec=[0 1 2 4]; restart = 0;
figure;
run('scattering_disk_penetrable',degree,PREC,nbEigVec,plotFlag,generalizedFlag,generalizedPMLFlag);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,PREC,nbEigVecList,plotFlag,generalizedFlag,generalizedPMLFlag)

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
    case {'none'}
        prec = 0;
    case {'CSL','ILU(CSL)','ILU(A)'}
        prec = 1;
    otherwise
        error('Error. \n%s is not a valid preconditioner', PREC);
end


[~, sysA] = computeSolNum2D_CG_heterogeneous(mesh, dofm, prec);


A = sysA.matA;
M = sysA.matM;
P = speye(size(A,2));
% b = sysA.rhsA;

% invPrec = @(x) P\x;
% switch PREC
%     case 'ILU(CSL)'
%         [L,U] = ilu(P);
%         invPrec = @(x) U\(L\x);
%     case 'ILU(A)'
%         [L,U] = ilu(A);
%         invPrec = @(x) U\(L\x);
% end


for iterDeflation = 1:size(nbEigVecList(:),1)
    nbEigVec = nbEigVecList(iterDeflation);


    % Compute deflation matrices
    if nbEigVec > 0
        switch benchmark
            case 'cavity'
                [eigvec,nbEigVec,Meigvec] = computeProjEigVec_cavity(mesh, dofm, nbEigVec,'closestEigvec',omega);
            case 'scattering_openCavity_NEU'
                [eigvec,nbEigVec,Meigvec] = computeProjEigVec_openCavity_NEU(mesh, dofm, nbEigVec, omega);
            case 'scattering_openCavity_DIR'
                [eigvec,nbEigVec,Meigvec] = computeProjEigVec_openCavity_DIR(mesh, dofm, nbEigVec, omega);
            case 'scattering_disk_penetrable'
                [eigvec,nbEigVec,Meigvec,eigvec_PML] = computeProjEigVec_het(mesh, dofm, nbEigVec);
            otherwise
                error('Error. \n%s is not a valid benchmark for deflation', benchmark);
        end
    [~,Pdef,~,~] = computeDefOp(nbEigVec, eigvec, A);

    if generalizedFlag
        [~,gPdef,~,~] = computeDefOp(nbEigVec, eigvec, M\A);
    end
    if generalizedPMLFlag
        [~,gPdefPML,~,~] = computeDefOp(nbEigVec, eigvec_PML, sysA.matMPML\A);
    end
    else
        Pdef = 1;
        if generalizedFlag
            gPdef = 1;
        end
        if generalizedPMLFlag
            gPdefPML = 1;
        end
    end

    [eigvec,eigval] = eigs(Pdef*A,10,'sm');
    eigval = diag(eigval);
    for i=1:length(eigval)
        filename = 'output/eigvec'+string(i)+'.pos';
        fieldname = 'eigvec'+string(i);
        writeField2D(dofm, mesh, eigvec(:,i), filename, fieldname);
    end
    
    if generalizedFlag
        [geigvec,geigval] = eigs(gPdef*M\A,10,'sm');
        geigval = diag(geigval);
        namefile = sprintf('output/geigval_%s_p%i_k=%g_prec=%s_def=%g.csv', benchmark, degree, omega, PREC, nbEigVec);
        writematrix(geigval, namefile, 'Delimiter', 'comma');
        for i=1:length(geigval)
            filename = 'output/geigvec'+string(i)+'.pos';
            fieldname = 'geigvec'+string(i);
            writeField2D(dofm, mesh, geigvec(:,i), filename, fieldname);
        end
    end
    if generalizedPMLFlag
        [geigvecPML,geigvalPML] = eigs(gPdefPML*sysA.matMPML\A,10,'sm');
        geigvalPML = diag(geigvalPML);
        namefile = sprintf('output/geigvalPML_%s_p%i_k=%g_prec=%s_def=%g.csv', benchmark, degree, omega, PREC, nbEigVec);
        writematrix(geigvalPML, namefile, 'Delimiter', 'comma');
        for i=1:length(geigvalPML)
            filename = 'output/geigvecPML'+string(i)+'.pos';
            fieldname = 'geigvecPML'+string(i);
            writeField2D(dofm, mesh, geigvecPML(:,i), filename, fieldname);
        end
    end

    namefile = sprintf('output/eigval_%s_p%i_k=%g_prec=%s_def=%g.csv', benchmark, degree, omega, PREC, nbEigVec);
    writematrix(eigval, namefile, 'Delimiter', 'comma');


    if plotFlag
        hold on;
        plot(real(eigval),imag(eigval),'o');
        if generalizedFlag
            plot(real(geigval),imag(geigval),'x');
        end
        xlabel('Re');
        ylabel('Im');
        title(sprintf('$\\sigma(\\mathbf{P}_{def}\\mathbf{A}), n_{def} = %g$', nbEigVec));
        legend('Location', 'best');
    end

end

% These variables must be cleared if 'cavity' is run after 'scattering'.
clear global LdomX LdomY LpmlX LpmlY;

end