%close all;
clear all;

degree = 3;
tau = 1;
BASIS = 1;
PREC = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE
benchmark = 'open';
h1 = log2(1/16); h2 = log2(1/34); hInt = (h2-h1)/3;
hList = 2.^[h1-2*hInt h1-hInt h1 h1+hInt h2-hInt h2 h2+hInt h2+2*hInt];
run(benchmark,degree,15*pi,hList,tau,BASIS,PREC);
run(benchmark,degree,30*pi,hList,tau,BASIS,PREC);

% BENCH CAVITY
benchmark = 'cavity';
h1 = log2(1/10); h2 = log2(1/15); hInt = (h2-h1)/2;
hList = 2.^[h1-hInt h1 h1+hInt h2 h2+hInt h2+2*hInt h2+3*hInt h2+4*hInt];
run(benchmark,degree,7.10*sqrt(2)*pi,hList,tau,BASIS,PREC);
run(benchmark,degree,7.01*sqrt(2)*pi,hList,tau,BASIS,PREC);

% BENCH WAVEGUIDE
benchmark = 'waveguide';
h1 = log2(1/8); h2 = log2(1/17); hInt = (h2-h1)/3;
hList = 2.^[h1-2*hInt h1-hInt h1 h1+hInt h2-hInt h2 h2+hInt h2+2*hInt];
run(benchmark,degree, 6*pi,hList,tau,BASIS,PREC);
run(benchmark,degree,12*pi,hList,tau,BASIS,PREC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,kList,hList,tau,BASIS,PREC)

disp(['---------------------------------------------------------']);
disp(['Method HDG - ' benchmark ' - k=' num2str(kList)]);
disp(['---------------------------------------------------------']);

global k;

Ndof = zeros(1,size(hList,2));
errorL2 = zeros(1,size(hList,2));
errorProjL2 = zeros(1,size(hList,2));
condGlo = zeros(1,size(hList,2));
condLocMin = zeros(1,size(hList,2));
condLocMax = zeros(1,size(hList,2));
invKH = zeros(1,size(hList,2));

for i = 1:size(hList,2)
    k = kList;
    h = hList(i);
    fprintf('%i/%i (h=%i)\n', i, size(hList,2), h)
    tic
    mesh = benchmark2D(benchmark,h);
    mesh = buildConnectivity2D(mesh);
    invKH(i) = 1/k*sqrt(mesh.numVer);  % mesh.numVer/4 for waveguide
    dofm = buildDofManager2D_DG(mesh, degree);
    Ndof(i) = dofm.numDofTRI;
    [solA, sysA, condLoc] = computeSolNum2D_HDG(mesh, dofm, tau, BASIS, PREC);
    errorL2(i) = computeNormError2D_DG(mesh, dofm, solA);
    solP = computeSolProjL2_2D_DG(mesh, dofm);
    errorProjL2(i) = computeNormError2D_DG(mesh, dofm, solP);
    fprintf('Errors:  %i  %i \n', errorL2(i), errorProjL2(i));
    condGlo(i) = condest(sysA.matS);
    fprintf('CondGlo:  %i \n', condGlo(i));
    condLocMin(i) = min(condLoc);
    condLocMax(i) = max(condLoc);
    fprintf('CondLoc:  %i  %i \n', condLocMin(i), condLocMax(i));
    toc
    disp('---------------------------------------------------------');
end

Dlambda = 2*pi/k * (sqrt(Ndof) - 1);

rezu1 = ["hList" "invKH" "Ndof" "Dlambda" "errorL2" "errorProjL2"];
rezu2 = [hList' invKH' Ndof' Dlambda' errorL2' errorProjL2'];
name = sprintf('output/errorVsH_HDG_%s_P%i_k%g_tau%g+%gi.csv', benchmark, degree, k, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

rezu1 = ["hList" "invKH" "Ndof" "Dlambda" "condGlo" "condLocMin" "condLocMax"];
rezu2 = [hList' invKH' Ndof' Dlambda' condGlo' condLocMin' condLocMax'];
name = sprintf('output/condVsH_HDG_%s_P%i_k%g_tau%g+%gi.csv', benchmark, degree, k, real(tau), imag(tau));
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

% figure;
% hold off;
% loglog(Dlambda, errorL2, '*-r');
% hold on;
% loglog(Dlambda, errorProjL2, '*:r');
% 
% figure;
% hold off;
% loglog(Dlambda, condGlo, '*-b');
% hold on;
% loglog(Dlambda, condLocMin, '*:r');
% loglog(Dlambda, condLocMax, '*-r');
% 
% figure;
% hold off;
% loglog(hList, condGlo, '*-b');
% hold on;
% loglog(hList, condLocMin, '*:r');
% loglog(hList, condLocMax, '*-r');

end