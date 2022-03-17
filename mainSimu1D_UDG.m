clear all;
close all;

addpath('benchmarks1D/');
addpath('kernels/');
addpath('tools/');

global k BCLeft BCRight

% Define parameters
degree = 3;
k = 40;
nume = 100;
numv = nume+1;
h = 1/nume;
tau = 1i; % 1i
resTol = 1e-4;
BCLeft = 'DIR';
BCRight = 'ABC';
TYPE = 'UDG-1'; % UDG-1 UDG-1b UDG-2

% Build mesh and dofManager
mesh = buildMesh1D(0, 1, numv);
dofm = buildDofManager1D_DG(mesh, degree);

disp(['=============================================']);
disp([TYPE ' - ' BCLeft '/' BCRight ' - degree ' num2str(degree) ' - k ' num2str(k) ' - nume ' num2str(nume) ' - tau ' num2str(tau)]);

% =========================================================================
% Compute solutions
% =========================================================================

% Compute numerical solution (full and reduced systems)
switch TYPE
    case 'UDG-1'
        [solFull, matA, rhsA, solRedu, matS, rhsS, matIIinv, matIG, rhsI] = ...
            computeSolNum1D_UDG1(mesh, dofm, @mySolP, @mySolU, @mySouP, @mySouU, tau);
    case 'UDG-1b'
        [solFull, matA, rhsA, solRedu, matS, rhsS, matIIinv, matIG, rhsI] = ...
            computeSolNum1D_UDG1b(mesh, dofm, @mySolP, @mySolU, @mySouP, @mySouU, tau);
    case 'UDG-2'
        [solFull, matA, rhsA, solRedu, matS, rhsS, matIIinv, matIG, rhsI] = ...
            computeSolNum1D_UDG2(mesh, dofm, @mySol, @mySolDer, @mySou, tau);
    otherwise
        warning('Error - No valid TYPE of equation.')
end


% Compute L2-projection solution
vecProjL2 = computeSolProjL2_1D_DG(mesh, dofm, @mySol);

% Vizu solution
%postProVizu1D_DG(mesh, dofm, @mySol, @mySolU, solRedu);

% =========================================================================
% Compute error
% =========================================================================

normL2sol       = computeNormSol1D_DG(mesh, @mySol);
normL2errSolNum = computeNormError1D_DG(mesh, dofm, @mySol, solRedu);
normL2errProjL2 = computeNormError1D_DG(mesh, dofm, @mySol, vecProjL2);
errSolNum = normL2errSolNum / normL2sol;
errProjL2 = normL2errProjL2 / normL2sol;
disp(['---------------------------------------------']);
disp(['    L2-Error           ' num2str(errSolNum)]);
disp(['    L2-Error (BestAp)  ' num2str(errProjL2)]);

% =========================================================================
% Iterative solution and spectrum
% =========================================================================

[vecSolIterA,~,~,iterA,resvecA] = gmres(matA,rhsA,size(matA,1),resTol,size(matA,1));
[vecSolIterS,~,~,iterS,resvecS] = gmres(matS,rhsS,size(matS,1),resTol,size(matS,1));
vecSolIterA = vecSolIterA(1:size(matA,2));
vecSolIterS = matIIinv*(rhsI-matIG*vecSolIterS);
normL2errSolIterA = computeNormError1D_DG(mesh, dofm, @mySol, vecSolIterA);
normL2errSolIterS = computeNormError1D_DG(mesh, dofm, @mySol, vecSolIterS);

[vecSolIterA,~,~,iterBiCGStabA,~] = bicgstab(matA,rhsA,resTol,10*size(matA,1));
[vecSolIterS,~,~,iterBiCGStabS,~] = bicgstab(matS,rhsS,resTol,10*size(matS,1));
vecSolIterA = vecSolIterA(1:size(matA,2));
vecSolIterS = matIIinv*(rhsI-matIG*vecSolIterS);
normL2errSolBiCGStabA = computeNormError1D_DG(mesh, dofm, @mySol, vecSolIterA);
normL2errSolBiCGStabS = computeNormError1D_DG(mesh, dofm, @mySol, vecSolIterS);

[vecSolIterA,~,~,iterCGNA,~] = conjgradn(matA,rhsA,resTol,size(matA,1));
[vecSolIterS,~,~,iterCGNS,~] = conjgradn(matS,rhsS,resTol,size(matS,1));
vecSolIterA = vecSolIterA(1:size(matA,2));
vecSolIterS = matIIinv*(rhsI-matIG*vecSolIterS);
normL2errSolCGNA = computeNormError1D_DG(mesh, dofm, @mySol, vecSolIterA);
normL2errSolCGNS = computeNormError1D_DG(mesh, dofm, @mySol, vecSolIterS);

[vecSolIterA,~,~,iterJacobiA,~] = jacobi(matA,rhsA,resTol,10*size(matA,1),0.5);
[vecSolIterS,~,~,iterJacobiS,~] = jacobi(matS,rhsS,resTol,10*size(matS,1),0.5);
vecSolIterA = vecSolIterA(1:size(matA,2));
vecSolIterS = matIIinv*(rhsI-matIG*vecSolIterS);
normL2errSolJacobiA = computeNormError1D_DG(mesh, dofm, @mySol, vecSolIterA);
normL2errSolJacobiS = computeNormError1D_DG(mesh, dofm, @mySol, vecSolIterS);

[eigenvecA,eigenvalA] = eigs(matA,size(matA,1));
[eigenvecS,eigenvalS] = eigs(matS,size(matS,1));
eigenvalA = diag(eigenvalA);
eigenvalS = diag(eigenvalS);
% singvalA = svds(matA,size(matA,1));
% singvalS = svds(matS,size(matS,1));

[eigenvecAA,eigenvalAA] = eigs(matA'*matA,size(matA,1));
[eigenvecSS,eigenvalSS] = eigs(matS'*matS,size(matS,1));
eigenvalAA = diag(eigenvalAA);
eigenvalSS = diag(eigenvalSS);

disp(['---------------------------------------------']);
disp(['A : Size               ' num2str(size(matA,1))]);
disp(['    Rank(eigenvectors) ' num2str(rank(eigenvecA))]);
disp(['    Cond(eigenvectors) ' num2str(cond(eigenvecA))]);
disp(['    Cond(A)            ' num2str(condest(matA))]);
disp(['    Cond(AA)           ' num2str(condest(matA'*matA))]);
disp(['    Det(A)             ' num2str(det(matA))]);
disp(['    Min eigenval AA    ' num2str(min(eigenvalAA))]);
disp(['    - - - - - - - - - - - - - - - - - - - - -']);
disp(['    IterGmres          ' num2str(iterA(2))]);
disp(['    Final L2-Error     ' num2str(normL2errSolIterA / normL2sol)]);
disp(['    IterBiCGS          ' num2str(iterBiCGStabA)]);
disp(['    Final L2-Error     ' num2str(normL2errSolBiCGStabA / normL2sol)]);
disp(['    IterCGN            ' num2str(iterCGNA)]);
disp(['    Final L2-Error     ' num2str(normL2errSolCGNA / normL2sol)]);
disp(['    IterRelax          ' num2str(iterJacobiA)]);
disp(['    Final L2-Error     ' num2str(normL2errSolJacobiA / normL2sol)]);
disp(['---------------------------------------------']);
disp(['S : Size               ' num2str(size(matS,1))]);
disp(['    Rank(eigenvectors) ' num2str(rank(eigenvecS))]);
disp(['    Cond(eigenvectors) ' num2str(cond(eigenvecS))]);
disp(['    Cond(S)            ' num2str(condest(matS))]);
disp(['    Cond(SS)           ' num2str(condest(matS'*matS))]);
disp(['    Det(S)             ' num2str(det(matS))]);
disp(['    Min eigenval SS    ' num2str(min(eigenvalSS))]);
disp(['    - - - - - - - - - - - - - - - - - - - - -']);
disp(['    IterGmres          ' num2str(iterS(2))]);
disp(['    Final L2-Error     ' num2str(normL2errSolIterS / normL2sol)]);
disp(['    IterBiCGS          ' num2str(iterBiCGStabS)]);
disp(['    Final L2-Error     ' num2str(normL2errSolBiCGStabS / normL2sol)]);
disp(['    IterCGN            ' num2str(iterCGNS)]);
disp(['    Final L2-Error     ' num2str(normL2errSolCGNS / normL2sol)]);
disp(['    IterRelax          ' num2str(iterJacobiS)]);
disp(['    Final L2-Error     ' num2str(normL2errSolJacobiS / normL2sol)]);
disp(['=============================================']);

disp(['\text{UDG} & full & ' ...
    num2str(errSolNum,'%.1e') ' & ' ...
    num2str(errProjL2,'%.1e') ' & ' ...
    num2str(size(matA,1)) ' & ' ...
    num2str(rank(eigenvecA)) ' & ' ...
    num2str(cond(eigenvecA),'%.1e') ' & ' ...
    num2str(condest(matA),'%.1e') ' & ' ...
    num2str(iterA(2)) ' & ' ...
    num2str(normL2errSolIterA / normL2sol,'%.1e') ' & ' ...
    num2str(iterBiCGStabA) ' & ' ...
    num2str(normL2errSolBiCGStabA / normL2sol,'%.1e') ' & ' ...
    num2str(iterCGNA) ' & ' ...
    num2str(normL2errSolCGNA / normL2sol,'%.1e') ' & ' ...
    num2str(iterJacobiA) ' & ' ...
    num2str(normL2errSolJacobiA / normL2sol,'%.1e') ' \\'
    ]);
disp(['\text{UDG} & red & ' ...
    num2str(errSolNum,'%.1e') ' & ' ...
    num2str(errProjL2,'%.1e') ' & ' ...
    num2str(size(matS,1)) ' & ' ...
    num2str(rank(eigenvecS)) ' & ' ...
    num2str(cond(eigenvecS),'%.1e') ' & ' ...
    num2str(condest(matS),'%.1e') ' & ' ...
    num2str(iterS(2)) ' & ' ...
    num2str(normL2errSolIterS / normL2sol,'%.1e') ' & ' ...
    num2str(iterBiCGStabS) ' & ' ...
    num2str(normL2errSolBiCGStabS / normL2sol,'%.1e') ' & ' ...
    num2str(iterCGNS) ' & ' ...
    num2str(normL2errSolCGNS / normL2sol,'%.1e') ' & ' ...
    num2str(iterJacobiS) ' & ' ...
    num2str(normL2errSolJacobiS / normL2sol,'%.1e') ' \\'
    ]);

% =========================================================================
% Vizu
% =========================================================================

% min(eigenvalAA) 2.4449e-05
% min(eigenvalSS) 2.4429e-04

% Eigenvalues and numerical range

% figure;
% hold off
% scatter(real(eigenvalA),imag(eigenvalA),'b','DisplayName','EigenVal A');
% hold on
% scatter(real(eigenvalS),imag(eigenvalS),'rx','DisplayName','EigenVal S');
% plot(fovals(matA,100),'-b','DisplayName','Numerical range A');
% plot(fovals(matS,100),'--r','DisplayName','Numerical range S');
% grid on; box on;
% title([TYPE ' (tau=' num2str(tau) ') - Eigenvalues : ' BCLeft ' + ' BCRight]);
% legend();
% 
% figure;
% hold off
% scatter(real(eigenvalAA),imag(eigenvalAA),'b','DisplayName','EigenVal AA');
% hold on
% scatter(real(eigenvalSS),imag(eigenvalSS),'rx','DisplayName','EigenVal SS');
% plot(fovals(matA'*matA,100),'-b','DisplayName','Numerical range AA');
% hold on
% plot(fovals(matS'*matS,100),'--r','DisplayName','Numerical range SS');
% grid on; box on;
% title([TYPE ' (tau=' num2str(tau) ') - Eigenvalues : ' BCLeft ' + ' BCRight]);
% legend();

% vecSolIterS = matIIinv*(rhsI-matIG*eigenvecS);
% for i=1:size(matS,1)
%     postProVizuFieldsDG(mesh, dofm, vecSolIterS(:,i));
%     title(eigenvalS(i));
% end

% figure;
% hold off
% semilogy(resvecA/resvecA(1),'-o','DisplayName','Residual UDG A');
% hold on
% semilogy(resvecS/resvecS(1),'-x','DisplayName','Residual UDG S');
% title('Residual history');
% legend();
