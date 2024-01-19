%close all;
clear all;

global k;
benchmark = 'open';
degree = 3;
k = 1;
hList = 2.^(-2:-0.5:-5);
tau = 1;

numVer = zeros(1,size(hList,2));
numTri = zeros(1,size(hList,2));
numEdg = zeros(1,size(hList,2));
numEdgBnd = zeros(1,size(hList,2));
numEdgInt = zeros(1,size(hList,2));
numDofPerLIN = zeros(1,size(hList,2));
numDofPerTRI = zeros(1,size(hList,2));

for i = 1:size(hList,2)
    h = hList(i);
    fprintf('   %i/%i (h=%i)\n', i, size(hList,2), hList(i));
    mesh = setupBenchmark2D(benchmark,h);
    mesh = buildMeshConnectivity(mesh);
    dofm = buildDofManager2D_DG(mesh, degree);
    
    numVer(i) = mesh.numVer;
    numTri(i) = mesh.numTri;
    numEdgBnd(i) = mesh.numEdgBnd;
    numEdgInt(i) = (3*mesh.numTri-mesh.numEdgBnd)/2;
    
    numDofPerLIN(i) = dofm.numDofPerLIN;
    numDofPerTRI(i) = dofm.numDofPerTRI;
    
    [~, sysDG] = computeSolNum2D_DG(mesh, dofm, tau, 1);
    [~, sysHDG] = computeSolNum2D_HDG(mesh, dofm, tau, 1);
    [~, sysCHDG] = computeSolNum2D_UDG(mesh, dofm, tau, 1);
    
    nDofDGnum(i) = size(sysDG.matA,1);
    nDofHDGnum(i) = size(sysHDG.matS,1);
    nDofCHDGnum(i) = size(sysCHDG.matS,1);
    nnzDGnum(i) = nnz(sysDG.matA);
    nnzHDGnum(i) = nnz(sysHDG.matS);
    nnzCHDGnum(i) = nnz(sysCHDG.matS);
    
end

numEdg = numEdgBnd+numEdgInt;

nDofDGtheo = 3*numTri.*numDofPerTRI;
nDofHDGtheo = numEdg.*numDofPerLIN;
nDofCHDGtheo = 3*numTri.*numDofPerLIN;
nnzDGtheo = numTri .* (7*numDofPerTRI.^2 + 54*numDofPerLIN.^2);
nnzHDGtheo = numEdg .* (5*numDofPerLIN.^2);
nnzCHDGtheo = 2*numEdg .* (3*numDofPerLIN.^2 + numDofPerLIN);

nnzDGtheoMin = numTri .* (7*numDofPerTRI.^2 + 54*numDofPerLIN.^2);
nnzHDGtheoMin = numEdg .* (5*numDofPerLIN.^2);
nnzCHDGtheoMin = 2*numEdg .* (3*numDofPerLIN.^2 + numDofPerLIN);

[nnzDGnum ; nnzDGtheo]
[nnzHDGnum ; nnzHDGtheo]
[nnzCHDGnum ; nnzCHDGtheo]

%Dlambda = 2*pi/hList * (sqrt(Ndof) - 1);

% rezu1 = ["khList" "condLoc"];
% rezu2 = [khList' condLocMaxHDG'];
% name = sprintf('output/condLocVsKH_HDG_%s_p%i_k%g.csv', benchmark, degree, k);
% writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');
% 
% rezu1 = ["khList" "condLoc"];
% rezu2 = [khList' condLocMaxUDG'];
% name = sprintf('output/condLocVsKH_CHDG_%s_p%i_k%g.csv', benchmark, degree, k);
% writematrix([rezu1 ; rezu2], name, 'Delimiter', 'semi');

figure(degree);
hold off;
loglog(1./hList, nnzDGtheo, 'b-');
hold on;
loglog(1./hList, nnzHDGtheo, 'r-');
loglog(1./hList, nnzCHDGtheo, 'g-');
loglog(1./hList, nnzDGnum, 'b*');
loglog(1./hList, nnzHDGnum, 'r*');
loglog(1./hList, nnzCHDGnum, 'g*');


