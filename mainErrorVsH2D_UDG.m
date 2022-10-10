%close all;
clear all;

headers2D;

tau = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE
benchmark = 'open'; degree = 3; kList = 15*pi; hList = 2.^(-(2:0.5:5));
run(benchmark,degree,kList,hList,tau);

% BENCH CAVITY
benchmark = 'cavity'; degree = 3; kList = (5+1/8)*sqrt(2)*pi; hList = 2.^(-(2:0.5:4));
run(benchmark,degree,kList,hList,tau);

% BENCH WAVEGUIDE
benchmark = 'waveguide'; degree = 3; kList = 6*pi; hList = 2.^(-(1:0.5:3.5)); tau = 1;
run(benchmark,degree,kList,hList,tau);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,kList,hList,tau)

global k;

Ndof = zeros(1,size(hList,2));
errorL2 = zeros(1,size(hList,2));
errorProjL2 = zeros(1,size(hList,2));
errorPostL2 = zeros(1,size(hList,2));
errorProjPostL2 = zeros(1,size(hList,2));
condGlo = zeros(1,size(hList,2));
condLocMin = zeros(1,size(hList,2));
condLocMax = zeros(1,size(hList,2));

for i = 1:size(hList,2)
    k = kList;
    h = hList(i);
    fprintf('   %i/%i (h=%i)\n', i, size(hList,2), h)
    tic
    mesh = benchmark2D(benchmark,h);
    mesh = buildMeshConnectivity(mesh);
    dofm = buildDofManager2D_DG(mesh, degree);
    Ndof(i) = dofm.numDofTRI;
    [solA, sysA, condLoc] = computeSolNum2D_UDG(mesh, dofm, tau, 1);
    [errorL2(i)] = computeNormError2D_DG(mesh, dofm, solA);
    [solP, ~] = computeSolProjL2_2D_DG(mesh, dofm);
    [errorProjL2(i)] = computeNormError2D_DG(mesh, dofm, solP);
    [solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solA);
    [errorPostL2(i)] = computeNormError2D_DG(mesh, dofmPost, solApost);
    [solPpost, ~] = computeSolProjL2_2D_DG(mesh, dofmPost);
    [errorProjPostL2(i)] = computeNormError2D_DG(mesh, dofmPost, solPpost);
    fprintf('   Errors:  %i  %i \n', errorL2(i), errorPostL2(i));
    condGlo(i) = condest(sysA.matS);
    condLocMin(i) = min(condLoc);
    condLocMax(i) = max(condLoc);
    fprintf('   CondGlo:  %i \n', condGlo(i));
    fprintf('   CondLoc:  %i  %i \n', condLocMin(i), condLocMax(i));
    toc
end

Dlambda = 2*pi/k * (sqrt(Ndof) - 1);

rezu1 = ["hList" "Ndof" "Dlambda" "errorL2" "errorProjL2" "errorPostL2" "errorProjPostL2"];
rezu2 = [hList' Ndof' Dlambda' errorL2' errorProjL2' errorPostL2' errorProjPostL2'];
name = sprintf('output/errorVsH_UDG_%s_P%i_k%g_tau%g+%gi.csv', benchmark, degree, k, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

rezu1 = ["hList" "Ndof" "Dlambda" "condGlo" "condLocMin" "condLocMax"];
rezu2 = [hList' Ndof' Dlambda' condGlo' condLocMin' condLocMax'];
name = sprintf('output/condVsH_UDG_%s_P%i_k%g_tau%g+%gi.csv', benchmark, degree, k, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

% figure(1);
% hold off;
% loglog(Dlambda, errorL2, '*-r');
% hold on;
% loglog(Dlambda, errorProjL2, '*:r');
% loglog(Dlambda, errorPostL2, '*-b');
% loglog(Dlambda, errorProjPostL2, '*:b');
% 
% figure(2);
% hold off;
% loglog(Dlambda, condGlo, '*-b');
% hold on;
% loglog(Dlambda, condLocMin, '*:r');
% loglog(Dlambda, condLocMax, '*-r');

end
