close all;
clear all;

global k R_disk L L_PML;

k = 25;
L = 1.1;
R_disk = 1;
L_PML = 0.1;

% benchmark = 'scatteringHard';
% benchmark = 'scatteringSoft';
benchmark = 'scatteringPML';


h = 0.01;
tol = 1e-10; maxit = 1000; itout = 50;

degree = 1;
PREC = 0;

mesh = benchmark2D(benchmark,h);
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


[solA, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);
errorL2 = computeNormError2D_CG(mesh, dofm, solA);
disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);

solP = computeSolProjL2_2D_CG(mesh, dofm);
errorL2 = computeNormError2D_CG(mesh, dofm, solP);
disp(['    L2-Error (refSol)   ' num2str(errorL2,'%1.2e')]);

writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
writeField2D(dofm, mesh, solP, 'output/solRef.pos', "solRef");
writeField2D(dofm, mesh, solA-solP, 'output/errNum.pos', "errNum");
system('gmsh output/solRef.pos output/solNum.pos output/errNum.pos&');

% [resVec, errorVec] = solverGMRES(mesh, dofm, sysA, tol, maxit, itout, @computeNormError2D_CG);


% sigma_x = zeros(dofm.numDofTRI,1);
% sigma_y = zeros(dofm.numDofTRI,1);

% for i=1:mesh.numVer
%     xQ = mesh.coord(i,1);
%     yQ = mesh.coord(i,2);
%     if abs(xQ) == 1 + h_PML
%         sigma_x(i) = 10^16;
%     elseif abs(yQ) == 1 + h_PML
%         sigma_y(i) = 10^16;
%     else
%         sigma_x(i) = 1./(1+h_PML-abs(xQ)) .* (1 <= abs(xQ)) .* (abs(xQ) <= 1 + h_PML);
%         sigma_y(i) = 1./(1+h_PML-abs(yQ)) .* (1 <= abs(yQ)) .* (abs(yQ) <= 1 + h_PML);
%     end
% end
% 
% writeField2D(dofm, mesh, sigma_x, 'output/sigma_x.pos', "sigma_x");
% writeField2D(dofm, mesh, sigma_y, 'output/sigma_y.pos', "sigma_y");
% system('gmsh output/sigma_x.pos output/sigma_y.pos&');







% X = linspace(-d,d,Npt);
% Y = linspace(-d,d,Npt);
% x = zeros(Npt,Npt);
% y = zeros(Npt,Npt);
% for i=1:Npt
%     x(i,:) = X(i);
%     y(:,i) = Y(i);
% end
% r = sqrt(x.^2+y.^2);
% 
% % p = solScattPlaneWaveHard(k,R,x,y);
% p = solScattPlaneWaveSoft(k,R,x,y);
% % 
% % mode = 20;
% % p = solScattModeHard(k,R,mode,x,y);
% % p = solScattModeSoft(k,R,mode,x,y);
% 
% p(r<R) = NaN;
% pcolor(x,y,real(p));
% shading interp;
% colorbar;