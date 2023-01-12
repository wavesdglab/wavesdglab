%close all;
clear all;

global k;
benchmark = 'open';
degree = 1;
h = 0.05;
tau = 1;

mesh = benchmark2D(benchmark,h);
mesh = buildMeshConnectivity(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

kList = 2.^(-8:1:0);
khList = kList*h;
khInvList = 1./khList;

condLocMaxHDG = zeros(1,size(kList,2));
condLocMaxUDG = zeros(1,size(kList,2));
for i = 1:size(kList,2)
    k = kList(i);
    fprintf('   %i/%i (k=%i)\n', i, size(khList,2), khList(i));
    [~, ~, condLocHDG] = computeSolNum2D_HDG(mesh, dofm, tau, 1);
    [~, ~, condLocUDG] = computeSolNum2D_UDG(mesh, dofm, tau, 1);
    condLocMaxHDG(i) = max(condLocHDG);
    condLocMaxUDG(i) = max(condLocUDG);
end

%Dlambda = 2*pi/kList * (sqrt(Ndof) - 1);

rezu1 = ["khList" "condLoc" "khListInv"];
rezu2 = [khList' condLocMaxHDG' khInvList'];
name = sprintf('output/condLocVsKH_HDG_%s_p%i_k%g.csv', benchmark, degree, k);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

rezu1 = ["khList" "condLoc" "khListInv"];
rezu2 = [khList' condLocMaxUDG' khInvList'];
name = sprintf('output/condLocVsKH_CHDG_%s_p%i_k%g.csv', benchmark, degree, k);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

% figure(degree);
% hold off;
% loglog(khList, condLocMaxHDG, 'r-*');
% hold on;
% loglog(khList, condLocMaxUDG, 'b-*');