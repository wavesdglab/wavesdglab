%close all;
clear all;

headers2D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE
p = 3;
k = 10*pi;
h = 1/16;
%C = h^(2*p)*k^(2*(p+1));
C = h^(2*p)*k^(2*p);
% benchmark = 'open'; degree = 3; kList = 2.^(2:0.5:6)*pi; hList = (C./kList.^(2*p)).^(1/(2*p)); tau = 1;
% run(benchmark,degree,C,kList,hList,tau);
benchmark = 'open'; degree = 3; kList = 2.^(2:0.5:6)*pi; hList = (C./kList.^(2*p)).^(1/(2*p)); tau = 1i;
run(benchmark,degree,C,kList,hList,tau);

% % BENCH CAVITY
% p = 3;
% k = (5+1/8)*sqrt(2)*pi;
% h = 1/8;
% %C = h^(2*p)*k^(2*(p+1));
% C = h^(2*p)*k^(2*p);
% benchmark = 'cavity'; degree = 3; kList = (4:0.25:7)*pi; hList = (C./kList.^(2*p)).^(1/(2*p)); tau = 1;
% run(benchmark,degree,C,kList,hList,tau);
% benchmark = 'cavity'; degree = 3; kList = (4:0.25:7)*pi; hList = (C./kList.^(2*p)).^(1/(2*p)); tau = 1i;
% run(benchmark,degree,C,kList,hList,tau);
% 
% % BENCH WAVEGUIDE
% p = 3;
% k = 6*pi;
% h = 1/8;
% %C = h^(2*p)*k^(2*(p+1));
% C = h^(2*p)*k^(2*p);
% benchmark = 'waveguide'; degree = 3; kList = 2.^(2:0.5:5)*pi; hList = (C./kList.^(2*p)).^(1/(2*p)); tau = 1;
% run(benchmark,degree,C,kList,hList,tau);
% benchmark = 'waveguide'; degree = 3; kList = 2.^(2:0.5:5)*pi; hList = (C./kList.^(2*p)).^(1/(2*p)); tau = 1i;
% run(benchmark,degree,C,kList,hList,tau);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,C,kList,hList,tau)

global k;

Ndof = zeros(1,size(kList,2));
errorL2 = zeros(1,size(kList,2));
errorProjL2 = zeros(1,size(kList,2));
errorPostL2 = zeros(1,size(kList,2));
errorProjPostL2 = zeros(1,size(kList,2));
condGlo = zeros(1,size(kList,2));
condLocMin = zeros(1,size(kList,2));
condLocMax = zeros(1,size(kList,2));

for i = 1:size(kList,2)
    k = kList(i);
    h = hList(i);
    fprintf('   %i/%i (k=%i)\n', i, size(kList,2), kList(i))
    tic
    mesh = benchmark2D(benchmark,h);
    mesh = buildMeshConnectivity(mesh);
    dofm = buildDofManager2D_DG(mesh, degree);
    Ndof(i) = dofm.numDofTRI;
    [solA, sysA, condLoc] = computeSolNum2D_UDG(mesh, dofm, tau);
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

rezu1 = ["kList" "hList" "Ndof" "Dlambda" "errorL2" "errorProjL2" "errorPostL2" "errorProjPostL2"];
rezu2 = [kList' hList' Ndof' Dlambda' errorL2' errorProjL2' errorPostL2' errorProjPostL2'];
name = sprintf('output/errorVsK_UDG_%s_P%i_C%g_tau%g+%gi.csv', benchmark, degree, C, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

rezu1 = ["kList" "hList" "Ndof" "Dlambda" "condGlo" "condLocMin" "condLocMax"];
rezu2 = [kList' hList' Ndof' Dlambda' condGlo' condLocMin' condLocMax'];
name = sprintf('output/condVsK_UDG_%s_P%i_C%g_tau%g+%gi.csv', benchmark, degree, C, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

% figure(1);
% hold off;
% loglog(kList, errorL2, '*-r');
% hold on;
% loglog(kList, errorProjL2, '*:r');
% loglog(kList, errorPostL2, '*-b');
% loglog(kList, errorProjPostL2, '*:b');
% 
% figure(2);
% hold off;
% loglog(kList, condGlo, '*-b');
% hold on;
% loglog(kList, condLocMin, '*:r');
% loglog(kList, condLocMax, '*-r');

end