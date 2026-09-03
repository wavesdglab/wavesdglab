close all;
clear;

global Options
Options.Basis = 'Jacobi';
Options.Error = 'Energy';

method = 'HDG';
degree = 3;
tau = 1;
BASIS = 1;
PREC = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BENCH FREE SPACE
benchmark = 'open';
h1 = log2(1/16); h2 = log2(1/34); hInt = (h2-h1)/3;
hList = 2.^[h1-2*hInt h1-hInt h1 h1+hInt h2-hInt h2 h2+hInt h2+2*hInt];
run(method,benchmark,degree,15*pi,hList,tau,BASIS,PREC);
run(method,benchmark,degree,30*pi,hList,tau,BASIS,PREC);

% BENCH CAVITY
benchmark = 'cavity';
h1 = log2(1/10); h2 = log2(1/15); hInt = (h2-h1)/2;
hList = 2.^[h1-hInt h1 h1+hInt h2 h2+hInt h2+2*hInt h2+3*hInt h2+4*hInt];
run(method,benchmark,degree,7.10*sqrt(2)*pi,hList,tau,BASIS,PREC);
run(method,benchmark,degree,7.01*sqrt(2)*pi,hList,tau,BASIS,PREC);

% BENCH WAVEGUIDE
benchmark = 'waveguide';
h1 = log2(1/8); h2 = log2(1/17); hInt = (h2-h1)/3;
hList = 2.^[h1-2*hInt h1-hInt h1 h1+hInt h2-hInt h2 h2+hInt h2+2*hInt];
run(method,benchmark,degree, 6*pi,hList,tau,BASIS,PREC);
run(method,benchmark,degree,12*pi,hList,tau,BASIS,PREC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(method,benchmark,degree,kList,hList,tau,BASIS,PREC)

disp(['---------------------------------------------------------']);
disp(['Method ' method ' - ' benchmark ' - k=' num2str(kList)]);
disp(['---------------------------------------------------------']);

global k h

hEff = zeros(1,size(hList,2));
invKH = zeros(1,size(hList,2));
NdofTRI = zeros(1,size(hList,2));
errorL2 = zeros(1,size(hList,2));
errorProjL2 = zeros(1,size(hList,2));
condGlo = zeros(1,size(hList,2));
specRad = zeros(1,size(hList,2));

for i = 1:size(hList,2)
    tic
    k = kList;
    h = hList(i);
    mesh = setupBenchmark2D(benchmark);
    mesh = buildConnectivity2D(mesh);
    hEff(i) = mesh.hmax;
    fprintf('%i/%i (h=%i ; heff=%i)\n', i, size(hList,2), h, hEff(i));
    invKH(i) = 1/(k*hEff(i));
    dofm = buildDofManager2D_DG(mesh, degree);
    NdofTRI(i) = dofm.numDofTRI;
    switch method
        case 'DG'
            theta = 1;
            [solA, sysA] = computeSolNum2D_DG(mesh, dofm, tau, theta, PREC);
        case 'HDG'
            [solA, sysA] = computeSolNum2D_HDG(mesh, dofm, tau, BASIS, PREC);
        case 'CHDG'
            [solA, sysA] = computeSolNum2D_CHDG(mesh, dofm, tau, BASIS, PREC);
            specRad(i) = max(abs(1-eigs(sysA.matS,100)));
            fprintf('1-SpecRad:   %i \n', 1-specRad(i));
    end
    errorL2(i) = computeNormError2D_DG(mesh, dofm, solA);
    fprintf('Error Num:   %i \n', errorL2(i));
    solP = computeSolProjL2_2D_DG(mesh, dofm);
    errorProjL2(i) = computeNormError2D_DG(mesh, dofm, solP);
    fprintf('Error Proj:  %i \n', errorProjL2(i));
    switch method
        case 'DG'
            condGlo(i) = condest(sysA.matA);
        case 'HDG'
            condGlo(i) = condest(sysA.matS);
        case 'CHDG'
            condGlo(i) = condest(sysA.matS);
    end
    fprintf('CondGlo:     %i \n', condGlo(i));
    toc
    disp('---------------------------------------------------------');
end

Dlambda = 2*pi/k * (sqrt(NdofTRI) - 1);

rezu1 = ["hList" "hEff" "invH" "invKH" "NdofTRI" "Dlambda" "errorL2" "errorProjL2" "condGlo"];
rezu2 = [hList' hEff' 1./hEff' invKH' NdofTRI' Dlambda' errorL2' errorProjL2' condGlo'];
name = sprintf('output/errorCondGloVsKH_%s_%s_P%i_k%g_tau%g+%gi_%g_%g.csv', method, benchmark, degree, k, real(tau), imag(tau), BASIS, PREC);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

% figure;
% loglog(Dlambda, errorL2, '*-r');
% hold on;
% loglog(Dlambda, errorProjL2, '*:r');
%
% figure;
% loglog(Dlambda, condGlo, '*-b');

end