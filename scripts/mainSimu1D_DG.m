clear all;
close all;

global k BCLeft BCRight

% Define parameters
degree = 3;
k = 40;
nume = 100;
numv = nume+1;
h = 1/nume;
theta = 1;
tau = 1; % 1i
resTol = 1e-4;
PREC = 'PrecNone'; % PrecNone PrecMass PrecMass2 PrecDiag
BCLeft = 'DIR';
BCRight = 'DIR';

% Build mesh and dofManager
mesh = buildMesh1D(0, 1, numv);
dofm = buildDofManager1D_DG(mesh, degree);

% -------------------------------------------------------------------------
% Compute solution and error
% -------------------------------------------------------------------------

disp(['---------------------------------------------------------']);
disp(['Method DG - Benchmark ' BCLeft '/' BCRight ]);
disp(['---------------------------------------------------------']);
disp(['    degree              ' num2str(degree)]);
disp(['    k                   ' num2str(k)]);
disp(['    nume                ' num2str(nume)]);
disp(['    theta               ' num2str(theta)]);
disp(['    tau                 ' num2str(tau)]);
disp(['    PREC                ' PREC]);
disp(['---------------------------------------------------------']);

[solA, sysA] = computeSolNum1D_DG(mesh, dofm, theta, tau, PREC);
errorL2 = computeNormError1D_DG(mesh, dofm, solA);

solP = computeSolProjL2_1D_DG(mesh, dofm);
errorProjL2 = computeNormError1D_DG(mesh, dofm, solP);

disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.6e')]);
disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.6e')]);
disp(['---------------------------------------------------------']);

% -------------------------------------------------------------------------
% Compute spectrum
% -------------------------------------------------------------------------

% [vecSolIterA,~,~,iterA,resvecA] = gmres(matA,rhsA,size(matA,1),resTol,size(matA,1));
% vecSolIterA = matP\vecSolIterA;
% normL2errSolIterA = computeNormError1D_DG(mesh, dofm, @mySolP, vecSolIterA);
% [eigenvecA,eigenvalA] = eigs(matA,size(matA,1));
% eigenvalA = diag(eigenvalA);
% singvalA = svds(matA,size(matA,1));
% 
% [vecSolIterA,~,~,iterBiCGStabA,~] = bicgstab(matA,rhsA,resTol,100*size(matA,1));
% vecSolIterA = matP\vecSolIterA;
% normL2errSolBiCGStabA = computeNormError1D_DG(mesh, dofm, @mySolP, vecSolIterA);
% 
% [vecSolIterA,~,~,iterCGNA,~] = conjgradn(matA,rhsA,resTol,size(matA,1));
% vecSolIterA = matP\vecSolIterA;
% normL2errSolCGNA = computeNormError1D_DG(mesh, dofm, @mySolP, vecSolIterA);
% 
% [vecSolIterA,~,~,iterJacobiA,~] = jacobi(matA,rhsA,resTol,10*size(matA,1),0.5);
% vecSolIterA = vecSolIterA(1:size(matA,2));
% normL2errSolJacobiA = computeNormError1D_DG(mesh, dofm, @mySolP, vecSolIterA);
% 
% disp(['---------------------------------------------']);
% disp(['A : Size               ' num2str(size(matA,1))]);
% disp(['    Eank(eigenvectors) ' num2str(rank(eigenvecA))]);
% disp(['    Cond(eigenvectors) ' num2str(cond(eigenvecA))]);
% disp(['    Cond(A)            ' num2str(condest(matA))]);
% disp(['    IterGmres          ' num2str(iterA(2))]);
% disp(['    Final L2-Error     ' num2str(normL2errSolIterA / normL2sol)]);
% disp(['    IterBiCGS          ' num2str(iterBiCGStabA)]);
% disp(['    Final L2-Error     ' num2str(normL2errSolBiCGStabA / normL2sol)]);
% disp(['    IterCGN            ' num2str(iterCGNA)]);
% disp(['    Final L2-Error     ' num2str(normL2errSolCGNA / normL2sol)]);
% disp(['    IterRelax          ' num2str(iterJacobiA)]);
% disp(['    Final L2-Error     ' num2str(normL2errSolJacobiA / normL2sol)]);
% disp(['=============================================']);
% 
% disp(['\text{DG} & & ' ...
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

% -------------------------------------------------------------------------
% Vizu
% -------------------------------------------------------------------------

% Eigenvalues and numerical range

% figure;
% hold off
% scatter(real(eigenvalA),imag(eigenvalA),'b','DisplayName','Eigenvalues');
% hold on
% plot(fovals(matA,100),'-b','DisplayName','Numerical range');
% grid on; box on;
% title([TYPE ' - Eigenvalues : ' BCLeft ' + ' BCRight]);
% legend();
%axis([-0.1 1.1 -1.5 1.2]);
%axis([-0.05 0.5 -1.5 1.2]);

% vecSolIterS = matIIinv*(rhsI-matIG*eigenvecS);
% for i=1:size(matS,1)
%     postProVizuFieldsDG(mesh, dofm, vecSolIterS(:,i));
%     title(eigenvalS(i));
% end

% figure(2);
% hold on
% semilogy(resvecA/resvecA(1),'-om','DisplayName','Residual DG A');
% title('Residual history');
% legend();

% figure(3);
% hold off
% scatter(real(singvalA),imag(singvalA),'DisplayName','SingVal DG A');
% grid on;
% title('Singular values');
% legend();

% figure(7);
% hold on
% scatter(real(eigenvalA),imag(eigenvalA),'DisplayName','EigenVal DG Upwind A');
% scatter(real(eigenvalA),imag(eigenvalA),'DisplayName','EigenVal DG Center A');
% grid on;
% title('Eigenvalues');
% axis([0 1 -1 1]);
% legend();

% figure(7);
% hold on
% scatter(real(eigenvalA),imag(eigenvalA),'DisplayName','ABC-ABC');
% grid on;
% title('Eigenvalues DG Upwind A');
% legend();


