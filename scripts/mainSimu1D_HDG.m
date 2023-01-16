clear all;
%close all;

global k BCLeft BCRight

% Define parameters
degree = 3;
k = 40;
numE = 200;
tau = 1; % 1i
resTol = 1e-4;
PREC = 'PrecNone'; % PrecNone PrecMass PrecMass2 PrecDiag
BCLeft = 'DIR';
BCRight = 'ABC';

% Build mesh and dofManager
mesh = buildMesh1D(0, 1, numE);
dofm = buildDofManager1D_DG(mesh, degree);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method HDG - Benchmark ' BCLeft '/' BCRight ]);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(1/numE)]);
disp(['    degree              ' num2str(degree)]);
disp(['    numE                ' num2str(numE)]);
disp(['    tau                 ' num2str(tau)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum1D_HDG(mesh, dofm, tau);
errorL2 = computeNormError1D_DG(mesh, dofm, solA);

solP = computeSolProjL2_1D_DG(mesh, dofm);
errorProjL2 = computeNormError1D_DG(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.6e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.6e')]);
disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Vizu solution
% -------------------------------------------------------------------------

plotField1D(mesh, dofm, solA, 'Numerical solution');

% -------------------------------------------------------------------------
% Iterative solution and spectrum
% -------------------------------------------------------------------------

% [vecSolIterA,~,~,iterA,resvecA] = gmres(matA,rhsA,size(matA,1),resTol,size(matA,1));
% [vecSolIterS,~,~,iterS,resvecS] = gmres(matS,rhsS,size(matS,1),resTol,size(matS,1));
% vecSolIterA = vecSolIterA(1:2*dofm.numDof);
% vecSolIterS = matIIinv*(rhsI-matIG*vecSolIterS);
% normL2errSolIterA = computeNormError1D_DG(mesh, dofm, @mySolP, vecSolIterA);
% normL2errSolIterS = computeNormError1D_DG(mesh, dofm, @mySolP, vecSolIterS);
% 
% [vecSolIterA,~,~,iterBiCGStabA,~] = bicgstab(matA,rhsA,resTol,100*size(matA,1));
% [vecSolIterS,~,~,iterBiCGStabS,~] = bicgstab(matS,rhsS,resTol,100*size(matS,1));
% vecSolIterA = vecSolIterA(1:2*dofm.numDof);
% vecSolIterS = matIIinv*(rhsI-matIG*vecSolIterS);
% normL2errSolBiCGStabA = computeNormError1D_DG(mesh, dofm, @mySolP, vecSolIterA);
% normL2errSolBiCGStabS = computeNormError1D_DG(mesh, dofm, @mySolP, vecSolIterS);
% 
% [vecSolIterA,~,~,iterCGNA,~] = conjgradn(matA,rhsA,resTol,100*size(matA,1));
% [vecSolIterS,~,~,iterCGNS,~] = conjgradn(matS,rhsS,resTol,100*size(matS,1));
% vecSolIterA = vecSolIterA(1:2*dofm.numDof);
% vecSolIterS = matIIinv*(rhsI-matIG*vecSolIterS);
% normL2errSolCGNA = computeNormError1D_DG(mesh, dofm, @mySolP, vecSolIterA);
% normL2errSolCGNS = computeNormError1D_DG(mesh, dofm, @mySolP, vecSolIterS);
% 
% [vecSolIterA,~,~,iterJacobiA,~] = jacobi(matA,rhsA,resTol,10*size(matA,1),0.5);
% [vecSolIterS,~,~,iterJacobiS,~] = jacobi(matS,rhsS,resTol,10*size(matS,1),0.5);
% vecSolIterA = vecSolIterA(1:size(matA,2));
% vecSolIterS = matIIinv*(rhsI-matIG*vecSolIterS);
% normL2errSolJacobiA = computeNormError1D_DG(mesh, dofm, @mySol, vecSolIterA);
% normL2errSolJacobiS = computeNormError1D_DG(mesh, dofm, @mySol, vecSolIterS);
% 
% [eigenvecA,eigenvalA] = eigs(matA,size(matA,1));
% [eigenvecS,eigenvalS] = eigs(matS,size(matS,1));
% eigenvalA = diag(eigenvalA);
% eigenvalS = diag(eigenvalS);
% % singvalA = svds(matA,size(matA,1));
% % singvalS = svds(matS,size(matS,1));
% 
% [eigenvecAA,eigenvalAA] = eigs(matA'*matA,size(matA,1));
% [eigenvecSS,eigenvalSS] = eigs(matS'*matS,size(matS,1));
% eigenvalAA = diag(eigenvalAA);
% eigenvalSS = diag(eigenvalSS);
% 
% disp(['---------------------------------------------']);
% disp(['A : Size               ' num2str(size(matA,1))]);
% disp(['    Rank(eigenvectors) ' num2str(rank(eigenvecA))]);
% disp(['    Cond(eigenvectors) ' num2str(cond(eigenvecA))]);
% disp(['    Cond(A)            ' num2str(condest(matA))]);
% disp(['    Cond(AA)           ' num2str(condest(matA'*matA))]);
% disp(['    Det(A)             ' num2str(det(matA))]);
% disp(['    Min eigenval AA    ' num2str(min(eigenvalAA))]);
% disp(['    - - - - - - - - - - - - - - - - - - - - -']);
% disp(['    IterGmres          ' num2str(iterA(2))]);
% disp(['    Final L2-Error     ' num2str(normL2errSolIterA / normL2sol)]);
% disp(['    IterBiCGS          ' num2str(iterBiCGStabA)]);
% disp(['    Final L2-Error     ' num2str(normL2errSolBiCGStabA / normL2sol)]);
% disp(['    IterCGN            ' num2str(iterCGNA)]);
% disp(['    Final L2-Error     ' num2str(normL2errSolCGNA / normL2sol)]);
% disp(['    IterRelax          ' num2str(iterJacobiA)]);
% disp(['    Final L2-Error     ' num2str(normL2errSolJacobiA / normL2sol)]);
% disp(['---------------------------------------------']);
% disp(['S : Size               ' num2str(size(matS,1))]);
% disp(['    Rank(eigenvectors) ' num2str(rank(eigenvecS))]);
% disp(['    Cond(eigenvectors) ' num2str(cond(eigenvecS))]);
% disp(['    Cond(S)            ' num2str(condest(matS))]);
% disp(['    Cond(SS)           ' num2str(condest(matS'*matS))]);
% disp(['    Det(S)             ' num2str(det(matS))]);
% disp(['    Min eigenval SS    ' num2str(min(eigenvalSS))]);
% disp(['    - - - - - - - - - - - - - - - - - - - - -']);
% disp(['    IterGmres          ' num2str(iterS(2))]);
% disp(['    Final L2-Error     ' num2str(normL2errSolIterS / normL2sol)]);
% disp(['    IterBiCGS          ' num2str(iterBiCGStabS)]);
% disp(['    Final L2-Error     ' num2str(normL2errSolBiCGStabS / normL2sol)]);
% disp(['    IterCGN            ' num2str(iterCGNS)]);
% disp(['    Final L2-Error     ' num2str(normL2errSolCGNS / normL2sol)]);
% disp(['    IterRelax          ' num2str(iterJacobiS)]);
% disp(['    Final L2-Error     ' num2str(normL2errSolJacobiS / normL2sol)]);
% disp(['=============================================']);
% 
% disp(['\text{HDG}(\tau=1) & full & ' ...
%     num2str(errSolNum,'%.1e') ' & ' ...
%     num2str(errProjL2,'%.1e') ' & ' ...
%     num2str(size(matA,1)) ' & ' ...
%     num2str(rank(eigenvecA)) ' & ' ...
%     num2str(cond(eigenvecA),'%.1e') ' & ' ...
%     num2str(condest(matA),'%.1e') ' & ' ...
%     num2str(iterA(2)) ' & ' ...
%     num2str(normL2errSolIterA / normL2sol,'%.1e') ' & ' ...
%     num2str(iterBiCGStabA) ' & ' ...
%     num2str(normL2errSolBiCGStabA / normL2sol,'%.1e') ' & ' ...
%     num2str(iterCGNA) ' & ' ...
%     num2str(normL2errSolCGNA / normL2sol,'%.1e') ' & ' ...
%     num2str(iterJacobiA) ' & ' ...
%     num2str(normL2errSolJacobiA / normL2sol,'%.1e') ' \\'
%     ]);
% disp(['\text{HDG}(\tau=1) & red & ' ...
%     num2str(errSolNum,'%.1e') ' & ' ...
%     num2str(errProjL2,'%.1e') ' & ' ...
%     num2str(size(matS,1)) ' & ' ...
%     num2str(rank(eigenvecS)) ' & ' ...
%     num2str(cond(eigenvecS),'%.1e') ' & ' ...
%     num2str(condest(matS),'%.1e') ' & ' ...
%     num2str(iterS(2)) ' & ' ...
%     num2str(normL2errSolIterS / normL2sol,'%.1e') ' & ' ...
%     num2str(iterBiCGStabS) ' & ' ...
%     num2str(normL2errSolBiCGStabS / normL2sol,'%.1e') ' & ' ...
%     num2str(iterCGNS) ' & ' ...
%     num2str(normL2errSolCGNS / normL2sol,'%.1e') ' & ' ...
%     num2str(iterJacobiS) ' & ' ...
%     num2str(normL2errSolJacobiS / normL2sol,'%.1e') ' \\'
%     ]);

% -------------------------------------------------------------------------
% Vizu
% -------------------------------------------------------------------------

% figure;
% hold off
% semilogy(resvecA/resvecA(1),'-o','DisplayName','Residual HDG A');
% hold on
% semilogy(resvecS/resvecS(1),'-x','DisplayName','Residual HDG S');
% title('Residual history');
% legend();

% figure;
% hold off
% scatter(real(singvalA),imag(singvalA),'DisplayName','SingVal HDG A');
% hold on
% scatter(real(singvalS),imag(singvalS),'x','DisplayName','SingVal HDG S');
% grid on;
% title('Singular values');
% legend();

% figure;
% hold off
% scatter(real(eigenvalA),imag(eigenvalA),'DisplayName','EigenVal HDG A');
% hold on
% scatter(real(eigenvalS),imag(eigenvalS),'x','DisplayName','EigenVal HDG S');
% grid on;
% title('Eigenvalues');
% legend();

% figure(5);
% hold off
% scatter(real(eigenvalA),imag(eigenvalA),'DisplayName','DIR-DIR');
% grid on;
% title('Eigenvalues HDG A');
% legend();
% 
% figure(6);
% hold off
% scatter(real(eigenvalS),imag(eigenvalS),'DisplayName','DIR-DIR');
% grid on;
% title('Eigenvalues HDG S');
% legend();