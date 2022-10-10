%close all;
clear all;

headers2D;

global k;
benchmark = 'open';
degree = 3;
kList = 2.^(-8:1:0);
h = 0.05;
khList = kList*h;
tau = 1;

Nver = zeros(1,size(kList,2));
Ndof = zeros(1,size(kList,2));
condLocMaxHDG = zeros(1,size(kList,2));
condLocMaxUDG = zeros(1,size(kList,2));

for i = 1:size(kList,2)
    k = kList(i);
    fprintf('   %i/%i (k=%i)\n', i, size(khList,2), khList(i));
    mesh = benchmark2D(benchmark,h);
    mesh = buildMeshConnectivity(mesh);
    dofm = buildDofManager2D_DG(mesh, degree);
    Ndof(i) = dofm.numDofTRI;
    Nver(i) = mesh.numVer;
    [~, ~, condLocHDG] = computeSolNum2D_HDG(mesh, dofm, tau, 1);
    [~, ~, condLocUDG] = computeSolNum2D_UDG(mesh, dofm, tau, 1);
    condLocMaxHDG(i) = max(condLocHDG);
    condLocMaxUDG(i) = max(condLocUDG);
end

%Dlambda = 2*pi/kList * (sqrt(Ndof) - 1);

rezu1 = ["khList" "condLoc"];
rezu2 = [khList' condLocMaxHDG'];
name = sprintf('output/condLocVsKH_HDG_%s_p%i_k%g.csv', benchmark, degree, k);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

rezu1 = ["khList" "condLoc"];
rezu2 = [khList' condLocMaxUDG'];
name = sprintf('output/condLocVsKH_CHDG_%s_p%i_k%g.csv', benchmark, degree, k);
writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

figure(degree);
hold off;
loglog(khList, condLocMaxHDG, 'r-*');
hold on;
loglog(khList, condLocMaxUDG, 'b-*');