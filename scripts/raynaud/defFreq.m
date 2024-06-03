
clear all;
%close all;

global k h

N=15;
LASTN = maxNumCompThreads(N);
disp(['---------------------------------------------------------']);
disp(['Previous maximum number of threads ' num2str(LASTN) ]);
disp(['Current maximum number of threads ' num2str(N) ]);
disp(['---------------------------------------------------------']);

computeSolNum2D = @computeSolNum2D_CG;

% Setup benchmark and parameters
benchmark = 'cavity';
switch benchmark
    case 'open'
        k = 15*pi;
        h = 1/16;
        tol = 1e-10; maxit = 1000; itout = 50;
    case 'cavity'
        k = 3.01*sqrt(2)*pi;
        h = 1/64;
        tol = 1e-6; maxit = 2000; itout =4;
        L = 1;
    case 'scatteringPML'
        k = 25;
        h = 0.05;
        tol = 1e-10; maxit = 2000; itout = 50;
        L = 1.1;
        R_disk = 1;
        L_PML = 0.2;
        computeSolNum2D = @computeSolNum2DPML_CG;
    case 'scattering_rec'
        global LdomX LdomY LpmlX LpmlY
        k = 10.5*pi;
        h = 1/8;
        tol = 1e-7; maxit = 5000; itout = 5;
        LdomX = 0.95;
        LdomY = 0.5;
        LpmlX = 0.8;
        LpmlY = 0.8;
    case 'waveguide'
        k = 6*pi;
        h = 1/8;
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 1; % P1
PREC = 0; % for preconditioner
% eigvecToDeflate = "closesteigvec"; %"firsteigvec" or "closesteigvec"
nbEigVec=1;

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree); % espace fonctionnel discret

tabfreq = 2.851*sqrt(2)*pi:0.001*sqrt(2)*pi:3.251*sqrt(2)*pi;
% rrGMRES = zeros(maxit, length(tabfreq));

% rrAD = zeros(maxit, length(tabfreq));
% rrD = zeros(maxit, length(tabfreq)); 
nbitGMRES = zeros(1, length(tabfreq));

nbitAD = zeros(1, length(tabfreq));
nbitD = zeros(1, length(tabfreq));

    [eigvec,nbEigVec] = computeEigVec2D_cavity_at_freq(mesh, dofm, nbEigVec,3.00001*sqrt(2)*pi);

for k = tabfreq
    
    
    Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1); % nb de points par longueur d'onde
    
    % -------------------------------------------------------------------------
    % Compute solution and error
    % -------------------------------------------------------------------------
    
    disp(['---------------------------------------------------------']);
    disp(['Method CG - Benchmark "' benchmark '"']);
    disp(['---------------------------------------------------------']);
    disp(['    k                   ' num2str(k)]);
    disp(['    h                   ' num2str(h)]);
    disp(['    degree              ' num2str(degree)]);
    disp(['    Dlambda             ' num2str(Dlambda)]);
    disp(['---------------------------------------------------------']);
    
    [~, sysA] = computeSolNum2D(mesh, dofm, PREC);
    
    A = sysA.matA;
    % M = sysA.matP;
    b = sysA.rhsA;
    
    if maxit > size(A,2)
        maxit = size(A,2);
    end

        [P,Q] = computeDefOp(nbEigVec, eigvec, A);
    
    %%%%%%%%%%% No deflation %%%%%%%%%%%
    
    % Compute GMRES without prec
    [xGMRES, ~, ~, itGMRES, rrGMRES] = gmres(A,b,[],tol,maxit);
    itGMRES = itGMRES(2);
    % rrGMRES = rrGMRES(:)./rrGMRES(1);
    % rrGMRES = rrGMRES(1:itout:end);
    
    % rrGMRES(1:length(rrGMRES), tabfreq == k) = rrGMRES;
    nbitGMRES(tabfreq == k) = itGMRES;
    
    
    

    

    
    
    
    
    
    % Compute GMRES with ADEF1 and closest eigvec : (P+Q)*A*x = (P+Q)*b
    [xAD, ~, ~, itAD, rrAD] = gmres((P+Q)*A,(P+Q)*b,[],tol,maxit);
    itAD = itAD(2);
    % rrAD = rrAD(:)./rrAD(1);
    % rrAD = rrAD(1:itout:end);
    
    % rrAD(1:length(rrAD), tabfreq == k) = rrAD;
    nbitAD(tabfreq == k) = itAD;
    
    
    % Compute GMRES with DEF1 and closest eigvec : P*A*x = P*b
    [xD, ~, ~, itD, rrD] = gmres(P*A,P*b,[],tol,maxit);
    itD = itD(2);
    % rrD = rrD(:)./rrD(1);
    % rrD = rrD(1:itout:end);
    
    nbitD(tabfreq == k) = itD;
    
    disp(['---------------------------------------------------------']);
    disp(['Number of iterations GMRES: ' num2str(itGMRES)]);
    disp(['Number of iterations GMRES with ADEF: ' num2str(nbitAD(tabfreq == k))]);
    disp(['Number of iterations GMRES with DEF: ' num2str(nbitD(tabfreq == k))]);
    disp(['---------------------------------------------------------']);
    
    
    
    
end

csvwrite('output/nbitGMRES.csv', nbitGMRES);
csvwrite('output/nbitAD.csv', nbitAD);
csvwrite('output/nbitD.csv', nbitD);

maxIt = max(max(nbitGMRES), max(max(nbitAD), max(nbitD)));
minIt = min(min(nbitGMRES), min(min(nbitAD), min(nbitD)));

figure
hold on
set(0,'DefaultFigureWindowStyle','docked')

p1 = semilogy(tabfreq,nbitGMRES,'b-o','DisplayName','rrGMRES','linewidth', 2,'markersize', 10);
p2 = semilogy(tabfreq,nbitAD,'r-+','DisplayName','rrAD','linewidth', 2,'markersize', 10);
p3 = semilogy(tabfreq,nbitD,'g-x','DisplayName','rrD','linewidth', 2,'markersize', 10);

set(gca, 'YScale', 'log')
box on
grid on
xlim([2.5*sqrt(2)*pi 3.5*sqrt(2)*pi]);
ylim auto;
title(['CG - ' benchmark ' - GMRES - k=' num2str(k) ' - h=' num2str(h) ' - degree=' num2str(degree) ' - nbEigvec=' num2str(nbEigVec)], 'interpreter', 'latex', 'fontsize', 20)
xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
legend('Location', 'southwest', 'fontsize', 15)


function [eigenvec,nbEigVec] = computeEigVec2D_cavity_at_freq(mesh, dofm, nbEigVec,freq)

mn = closest_to_k(nbEigVec, freq);

% Quadrature and shape functions
degreeQ = 4*dofm.degree;
[uQ, vQ, weights] = quadratureGaussTRI(degreeQ);
weights = sparse(1:size(weights,1), 1:size(weights,1), weights);
shapeQ = functionsShapeTRI(uQ, vQ, dofm.degree);

% Build matrix and RHS vector
matP = sparse(dofm.numDofTRI, dofm.numDofTRI);
rhsP = zeros(dofm.numDofTRI, size(mn, 1));

for tri=1:mesh.numTri
    
    % Mapping
    ver = mesh.mapTriToVer(tri,:);
    V1 = mesh.coord(ver(1),:);
    V2 = mesh.coord(ver(2),:);
    V3 = mesh.coord(ver(3),:);
    [xQ, yQ] = locToGloTRI(uQ, vQ, V1, V2, V3);
    Jdxdu = [(V2-V1)' (V3-V1)'] * 0.5;  % [ dx/du dx/dv ; dy/du dy/dv ]
    detJdxdu = abs(det(Jdxdu));
    
    % Reference solution
    refQ = zeros(size(xQ,1), size(mn, 1));
    for i=1:size(mn, 1)
        m = mn(i, 1);
        n = mn(i, 2);
        
        refQ(:, i) = (sin(m*xQ*pi).*sin(n*yQ*pi)) .* (16*m*n*pi^2 /(pi^2*(m^2 + n^2)-freq^2));
    end
    
    % Orientation
    orientation = ones(dofm.numDofPerTRI,1);
    if(ver(1) > ver(2))
        orientation(dofm.locEdg(1,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(ver(2) > ver(3))
        orientation(dofm.locEdg(2,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    if(ver(3) > ver(1))
        orientation(dofm.locEdg(3,:)) = (-1).^(0:dofm.numDofPerEdg-1);
    end
    orientation = sparse(1:dofm.numDofPerTRI, 1:dofm.numDofPerTRI, orientation);
    
    % Shape functions with orientation
    shapeOrQ = shapeQ * orientation;
    
    % Local matrix and RHS vector
    matPel = shapeOrQ' * weights * shapeOrQ * detJdxdu;
    rhsPel = shapeOrQ' * weights * refQ * detJdxdu;
    
    % Assembling
    dof = dofm.locToGloTRI(tri,:);
    matP(dof,dof) = matP(dof,dof) + matPel;
    rhsP(dof,:) = rhsP(dof,:) + rhsPel;
    
    
    nbEigVec = size(mn, 1);
    
end

% Solution
eigenvec = matP\rhsP;

% for i=1:size(mn, 1)
%     filename = 'output/eigenvec' + string(i) + '.pos';
%     fieldname = 'eigenvec' + string(i);
%     writeField2D(dofm, mesh, eigenvec(:,i), filename, fieldname);
%     fprintf('field [%i] saved \n', i);
% end

end

function indices = smallest_sum_of_squares(nb)
% Generate a grid of integer indices up to a certain limit
limit = 100;
[M, N] = meshgrid(1:limit, 1:limit);

% Compute the sum of squares m^2 + n^2
sum_squares = M.^2 + N.^2;

% Flatten the matrices and sort the indices based on sum of squares
[~, sorted_indices] = sort(sum_squares(:));

% Extract the nb+1 smallest indices
nb_smallest_indices = sorted_indices(1:nb+1);

% Convert linear indices to subscripts
[m, n] = ind2sub(size(sum_squares), nb_smallest_indices);

%% test if the nbth and the nb+1th smallest are equal
if sum_squares(m(nb), n(nb)) == sum_squares(m(nb+1), n(nb+1))
    indices = [m(1:nb+1), n(1:nb+1)];
else
    indices = [m(1:nb), n(1:nb)];
end
end

function indices = closest_to_k(nb, k)
% Generate a grid of integer indices up to a certain limit
limit = 100;
[M, N] = meshgrid(1:limit, 1:limit);

% Compute the difference between m^2 + n^2 and k^2/pi^2
diff = abs(M.^2 + N.^2 - k^2/pi^2);

% Flatten the matrices and sort the indices based on the difference
[~, sorted_indices] = sort(diff(:));

% Extract the nb+1 smallest indices
indices = sorted_indices(1:nb+1);

% Convert linear indices to subscripts
[m, n] = ind2sub(size(diff), indices);

%% test if the nbth and the nb+1th smallest are equal
if diff(m(nb), n(nb)) == diff(m(nb+1), n(nb+1))
    indices = [m(1:nb+1), n(1:nb+1)];
else
    indices = [m(1:nb), n(1:nb)];
end
end