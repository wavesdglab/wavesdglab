close all;
clear all;

global k h

benchmark = 'open';
h = 0.04;
tau = 1;
BASIS = 1;
PREC = 0;

mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
h = mesh.hmax;

kList = 2.^(-8:1:0);
khList = kList*h;
khInvList = 1./khList;

figure;
for degree = 1:3
    
    dofm = buildDofManager2D_DG(mesh, degree);
    
    condLocMaxHDG = zeros(1,size(kList,2));
    condLocMaxCHDG = zeros(1,size(kList,2));
    for i = 1:size(kList,2)
        k = kList(i);
        fprintf('   %i/%i (k=%i)\n', i, size(khList,2), khList(i));
        [~, ~, condLocHDG] = computeSolNum2D_HDG(mesh, dofm, tau, 1, 0);
        [~, ~, condLocCHDG] = computeSolNum2D_CHDG(mesh, dofm, tau, 1, 0);
        condLocMaxHDG(i) = max(condLocHDG);
        condLocMaxCHDG(i) = max(condLocCHDG);
    end
    
    rezu1 = ["khList" "condLoc" "khListInv"];
    rezu2 = [khList' condLocMaxHDG' khInvList'];
    name = sprintf('output/condLocVsKH_HDG_p%i_k%g_tau%g+%gi.csv', degree, k, real(tau), imag(tau));
    writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');
    
    rezu1 = ["khList" "condLoc" "khListInv"];
    rezu2 = [khList' condLocMaxCHDG' khInvList'];
    name = sprintf('output/condLocVsKH_CHDG_p%i_k%g_tau%g+%gi.csv', degree, k, real(tau), imag(tau));
    writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');
    
    loglog(khList, condLocMaxHDG, 'r-*');
    hold on;
    loglog(khList, condLocMaxCHDG, 'b-*');
end