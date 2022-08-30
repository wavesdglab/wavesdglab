%close all;
clear all;

headers2D;
global k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE [ P1 10pi h=1/32 ; P3 40pi h=1/32 ]
% k = 20*pi;
% hList = 2.^(-(1:0.5:6));
% degree = 3;
% benchmark = 'open';

% BENCH WAVEGUIDE [ P1 2pi h=1/16 ; P3 6pi h=1/8 ]
% k = 6*pi;
% hList = 2.^(-(0.5:0.5:4.5));
% degree = 3;
% benchmark = 'waveguide';

% BENCH CAVITY [ P1 (3*sqrt(2)*pi+sqrt(2)*pi/8 or /64) h=1/32 ; P3 (5*sqrt(2)*pi+sqrt(2)*pi/8 or /64) h=1/8 ]
k = 5*sqrt(2)*pi+sqrt(2)*pi/64;
hList = 2.^(-(1:0.5:5));
degree = 3;
benchmark = 'cavity';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

rezu = [hList' Dlambda' errorL2' errorProjL2' errorPostL2' errorProjPostL2'];
writematrix(rezu,'output/errorVsH_UDG_Open_1.csv','Delimiter','semi');

figure(1);
hold off;
loglog(Dlambda, errorL2, '*-r');
hold on;
loglog(Dlambda, errorProjL2, '*:r');
loglog(Dlambda, errorPostL2, '*-b');
loglog(Dlambda, errorProjPostL2, '*:b');
axis([1 100 1e-8 10])

% rezu = ones(5);
% writematrix(rezu,'test.csv','Delimiter','tab')
