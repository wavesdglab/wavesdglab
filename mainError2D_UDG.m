%close all;
clear all;

headers2D;
global k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE [ P1 10pi h=1/32 ; P3 40pi h=1/32 ]
% benchmark = 'open'; degree = 1; k = 10*pi; hList = 2.^(-(2.5:0.5:6));
% benchmark = 'open'; degree = 3; k = 40*pi; hList = 2.^(-(2.5:0.5:6));

% BENCH WAVEGUIDE [ P1 2pi h=1/16 ; P3 6pi h=1/8 ]
% benchmark = 'waveguide'; degree = 1; k = 2*pi; hList = 2.^(-(1:0.5:5));
benchmark = 'waveguide'; degree = 3; k = 6*pi; hList = 2.^(-(0.5:0.5:4.5));

% BENCH CAVITY [ P1 (3*sqrt(2)*pi+sqrt(2)*pi/8 or /64) h=1/32 ; P3 (5*sqrt(2)*pi+sqrt(2)*pi/8 or /64) h=1/8 ]
% benchmark = 'cavity'; degree = 1; k = 3*sqrt(2)*pi+sqrt(2)*pi/8; hList = 2.^(-(2:0.5:5));
% benchmark = 'cavity'; degree = 3; k = 5*sqrt(2)*pi+sqrt(2)*pi/8; hList = 2.^(-(1:0.5:4));

tau = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndof = zeros(1,size(hList,2));
errorL2 = zeros(1,size(hList,2));
errorProjL2 = zeros(1,size(hList,2));
errorPostL2 = zeros(1,size(hList,2));
errorProjPostL2 = zeros(1,size(hList,2));

for i = 1:size(hList,2)
    fprintf('   %i/%i (h=%i)\n', i, size(hList,2), hList(i))
    tic
    mesh = benchmark2D(benchmark,hList(i));
    mesh = buildMeshConnectivity(mesh);
    dofm = buildDofManager2D_DG(mesh, degree);
    Ndof(i) = dofm.numDofTRI;
    [solA, sysA] = computeSolNum2D_UDG(mesh, dofm, tau);
    [errorL2(i)] = computeNormError2D_DG(mesh, dofm, solA);
    [solP, sysP] = computeSolProjL2_2D_DG(mesh, dofm);
    [errorProjL2(i)] = computeNormError2D_DG(mesh, dofm, solP);
    [solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solA);
    [errorPostL2(i)] = computeNormError2D_DG(mesh, dofmPost, solApost);
    [solPpost, sysPpost] = computeSolProjL2_2D_DG(mesh, dofmPost);
    [errorProjPostL2(i)] = computeNormError2D_DG(mesh, dofmPost, solPpost);
    fprintf('   Errors:  %i  %i \n', errorL2(i), errorPostL2(i));
    toc
end

Dlambda = 2*pi/k * (sqrt(Ndof) - 1);

rezu1 = ["hList" "Ndof" "Dlambda" "errorL2" "errorProjL2" "errorPostL2" "errorProjPostL2"];
rezu2 = [hList' Ndof' Dlambda' errorL2' errorProjL2' errorPostL2' errorProjPostL2'];
name = sprintf('output/errorVsH_UDG_%s_P%i_k%g_tau%g+%gi.csv', benchmark, degree, k, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

figure(1);
hold off;
loglog(Dlambda, errorL2, '*-r');
hold on;
loglog(Dlambda, errorProjL2, '*:r');
loglog(Dlambda, errorPostL2, '*-b');
loglog(Dlambda, errorProjPostL2, '*:b');