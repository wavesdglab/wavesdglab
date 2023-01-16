clear all;
close all;

global k BCLeft BCRight

% Setup benchmark and parameters
degree = 3;
k = 40;
nume = 200;
numv = nume+1;
h = 1/nume;
resTol = 1e-4;
PREC = 'PrecNone'; % PrecNone PrecMass PrecDiag PrecShiftLap
alphaPrec = 1; % 1+1i
BCLeft = 'DIR';
BCRight = 'ABC';

% Build mesh and DOF manager
mesh = buildMesh1D(0, 1, numv);
dofm = buildDofManager1D_CG(mesh, degree);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method CG - Benchmark ' BCLeft '/' BCRight ]);
disp(['---------------------------------------------------------']);
disp(['    k                   ' num2str(k)]);
disp(['    h                   ' num2str(h)]);
disp(['    degree              ' num2str(degree)]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum1D_CG(mesh, dofm, PREC, alphaPrec);
errorL2 = computeNormError1D_CG(mesh, dofm, solA);

solP = computeSolProjL2_1D_CG(mesh, dofm);
errorProjL2 = computeNormError1D_CG(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.6e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.6e')]);
disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Vizu solution
% -------------------------------------------------------------------------

% plotSol1D_CG(mesh, dofm, solA, 'Numerical solution');

% -------------------------------------------------------------------------
% Iterative solution and spectrum
% -------------------------------------------------------------------------

% [matM, matK, matD] = buildMatrixGlo1D_CG(mesh, dofm);
% [vecSolIterA,~,~,iterA,resvecA] = gmres(matA,rhsA,size(matA,1),resTol,size(matA,1));
% [vecSolIterS,~,~,iterS,resvecS] = gmres(matS,rhsS,size(matS,1),resTol,size(matS,1));
% solG = vecSolIterS;
% solI = matII\(rhsI-matIG*solG);
% vecSolIterS = [ solG ; solI ];
% vecSolIterA = matP\vecSolIterA;
% vecSolIterS = matP\vecSolIterS;
% [normL2errSolIterA,~,~] = computeNormError1D_CG(mesh, dofm, @mySol, @mySolDer, vecSolIterA);
% [normL2errSolIterS,~,~] = computeNormError1D_CG(mesh, dofm, @mySol, @mySolDer, vecSolIterS);
% 
% [vecSolIterA,~,~,iterCGNA,~] = solverCGN(matA,rhsA,resTol,size(matA,1));
% [vecSolIterS,~,~,iterCGNS,~] = solverCGN(matS,rhsS,resTol,size(matS,1));
% solG = vecSolIterS;
% solI = matII\(rhsI-matIG*solG);
% vecSolIterS = [ solG ; solI ];
% vecSolIterA = matP\vecSolIterA;
% vecSolIterS = matP\vecSolIterS;
% [normL2errSolCGNA,~,~] = computeNormError1D_CG(mesh, dofm, @mySol, @mySolDer, vecSolIterA);
% [normL2errSolCGNS,~,~] = computeNormError1D_CG(mesh, dofm, @mySol, @mySolDer, vecSolIterS);
% 
% [vecSolIterA,~,~,iterJacobiA,~] = jacobi(matA,rhsA,resTol,10*size(matA,1),0.5);
% [vecSolIterS,~,~,iterJacobiS,~] = jacobi(matS,rhsS,resTol,10*size(matS,1),0.5);
% solG = vecSolIterS;
% solI = matII\(rhsI-matIG*solG);
% vecSolIterS = [ solG ; solI ];
% vecSolIterA = matP\vecSolIterA;
% vecSolIterS = matP\vecSolIterS;
% [normL2errSolJacobiA,~,~] = computeNormError1D_CG(mesh, dofm, @mySol, @mySolDer, vecSolIterA);
% [normL2errSolJacobiS,~,~] = computeNormError1D_CG(mesh, dofm, @mySol, @mySolDer, vecSolIterS);
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
% % matM = buildMatrixGloMassCG(mesh, dofm);
% % matK = buildMatrixGloStiffnessCG(mesh, dofm);
% % [eigenvecM,eigenvalM] = eigs(matM,size(matM,1));
% % [eigenvecK,eigenvalK] = eigs(matK,size(matmatKM,1));
% % eigenvalM = diag(eigenvalM);
% % eigenvalK = diag(eigenvalK);
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
% disp(['    IterCGN            ' num2str(iterCGNS)]);
% disp(['    Final L2-Error     ' num2str(normL2errSolCGNS / normL2sol)]);
% disp(['    IterRelax          ' num2str(iterJacobiS)]);
% disp(['    Final L2-Error     ' num2str(normL2errSolJacobiS / normL2sol)]);
% disp(['=============================================']);
% 
% disp(['\text{CG} & full & ' ...
%     num2str(errSol,'%.1e') ' & ' ...
%     num2str(errSolProjL2,'%.1e') ' & ' ...
%     num2str(size(matA,1)) ' & ' ...
%     num2str(rank(eigenvecA)) ' & ' ...
%     num2str(cond(eigenvecA),'%.1e') ' & ' ...
%     num2str(condest(matA),'%.1e') ' & ' ...
%     num2str(iterA(2)) ' & ' ...
%     num2str(normL2errSolIterA / normL2sol,'%.1e') ' & ' ...
%     num2str(iterCGNA) ' & ' ...
%     num2str(normL2errSolCGNA / normL2sol,'%.1e') ' & ' ...
%     num2str(iterJacobiA) ' & ' ...
%     num2str(normL2errSolJacobiA / normL2sol,'%.1e') ' \\'
%     ]);
% disp(['\text{CG} & red & ' ...
%     num2str(errSol,'%.1e') ' & ' ...
%     num2str(errSolProjL2,'%.1e') ' & ' ...
%     num2str(size(matS,1)) ' & ' ...
%     num2str(rank(eigenvecS)) ' & ' ...
%     num2str(cond(eigenvecS),'%.1e') ' & ' ...
%     num2str(condest(matS),'%.1e') ' & ' ...
%     num2str(iterS(2)) ' & ' ...
%     num2str(normL2errSolIterS / normL2sol,'%.1e') ' & ' ...
%     num2str(iterCGNS) ' & ' ...
%     num2str(normL2errSolCGNS / normL2sol,'%.1e') ' & ' ...
%     num2str(iterJacobiS) ' & ' ...
%     num2str(normL2errSolJacobiS / normL2sol,'%.1e') ' \\'
%     ]);

% -------------------------------------------------------------------------
% Vizu
% -------------------------------------------------------------------------

% Eigenvalues and numerical range

% figure;
% hold off
% scatter(real(eigenvalA),imag(eigenvalA),'b','DisplayName','EigenVal A');
% hold on
% scatter(real(eigenvalS),imag(eigenvalS),'rx','DisplayName','EigenVal S');
% plot(fovals(matA,300),'-b','DisplayName','Numerical range A');
% plot(fovals(matS,300),'--r','DisplayName','Numerical range S');
% grid on; box on;
% title(['CG' ' - Eigenvalues : ' BCLeft ' + ' BCRight]);
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
% title(['CG' ' - Eigenvalues : ' BCLeft ' + ' BCRight]);
% legend();


% figure(7)
% hold on
% semilogy(resvecA/resvecA(1),'-o','DisplayName','Residual CG A');

% figure(6)
% hold on
% semilogy(resvecS/resvecS(1),'-o','DisplayName','Residual CG S');

% figure(2);
% hold off
% semilogy(resvecA/resvecA(1),'-o','DisplayName','Residual CG A');
% hold on
% semilogy(resvecS/resvecS(1),'-x','DisplayName','Residual CG S');
% title('Residual history');
% legend();
%
% figure(3);
% hold off
% scatter(real(singvalA),imag(singvalA),'DisplayName','SingVal CG A');
% hold on
% scatter(real(singvalS),imag(singvalS),'x','DisplayName','SingVal CG S');
% grid on;
% title('Singular values');
% legend();

% figure(1);
% hold on;
% scatter(real(eigenvalA),imag(eigenvalA),'DisplayName','DIR-DIR');
% grid on;
% box on;
% title('Eigenvalues CG A');
% legend();

% Z = -1000:0.01:0;
% circ = (1+Z)./(alphaPrec+Z);
% circX = real(circ);
% circY = imag(circ);
% hold on
% plot(circX,circY)

% figure(2);
% [eigenval,listOrder] = sort(real(eigenvalA));
% for i=listOrder'
%     plotSolCG(mesh, dofm, eigenvecA(:,i), ['Eigenvec ' num2str(i) ' : ' num2str(eigenvalA(i))]);
%     pause;
% end

% figure(2);
% hold on;
% scatter(real(eigenvalS),imag(eigenvalS),'DisplayName','DIR-DIR');
% grid on;
% box on;
% title('Eigenvalues CG S');
% legend();
