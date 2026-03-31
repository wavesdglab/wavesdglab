%% Setup
%close all;
clear;
clear global;

global k h Options

Options.Basis = 'Jacobi'; % Jacobi, Lagrange
Options.Error = 'L2'; % L2, H1


%%%%%%%%%%%%%%%%%%%%%%%% Cavity

k = 3.01*sqrt(2)*pi;
h = 1./20;
tol = 1e-10;
degree = 1;

mesh = setupBenchmark2D('cavity'); % cavity or scattering_openCavity_NEU
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);
Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

[~, sysA] = computeSolNum2D_CG_red(mesh, dofm, 'none');
A = sysA.matA;
P = sysA.matP;
M = sysA.matM;
Melim = sysA.matMelim;
M0 = sysA.matMevpb;
M0_inv = M0;
nonZeroRows = any(M0, 2);
M0_inv(nonZeroRows, nonZeroRows) = inv(M0(nonZeroRows, nonZeroRows));
I_int = sysA.matI_int;
b = sysA.rhsA;
I = eye(size(A,2));
R_int = sysA.matR_int;


Ared = sysA.matAred;
Mred = sysA.matMred;
bred = sysA.rhsAred;

nbEigVec = 20;
[u,nbEigVec,w,matP] = computeProjEigVec_cavity_bis(mesh, dofm, nbEigVec,'closestEigvec',k);
% [u,nbEigVec] = computeProjEigVec_openCavity_NEU(mesh, dofm, nbEigVec, k);
[v,~] = eigs(A,nbEigVec,'sm');
% w = Melim*u;
% Afun = @(x) (A/M)'*x;
% [vp,~] = eigs(Afun,nbEigVec,'sm'); % ((AM^(-1))'*vp = lambda*vp) <=> (vp'*A = lambda vp'*M)

uu = u;
vv = v;
command = strcat('gmsh output/mesh.msh');
for i = 1:nbEigVec
    L2norm = sqrt(u(:,i)'*Melim*u(:,i));
    u(:,i) = u(:,i)./L2norm;
    l2norm = sqrt(uu(:,i)'*uu(:,i));
    uu(:,i) = uu(:,i)./l2norm;
%     u(:,i) = u(:,i)/max(u(:,i));
%     namefile = sprintf('output/projvec_%g.pos',i);
%     namesol = strcat('projvec_',num2str(i));
%     writeField2D(dofm, mesh, u(:,i), namefile, namesol);

    L2norm = sqrt(v(:,i)'*Melim*v(:,i));
    v(:,i) = v(:,i)./L2norm;
    l2norm = sqrt(vv(:,i)'*vv(:,i));
    vv(:,i) = vv(:,i)./l2norm;
% %     v(:,i) = v(:,i)/max(v(:,i));
%     namefile = sprintf('output/eigvec_%g.pos',i);
%     namesol = strcat('eigvec_',num2str(i));
%     writeField2D(dofm, mesh, w(:,i), namefile, namesol);
%     command = strcat(command, ' output/projvec_',num2str(i),'.pos output/eigvec_',num2str(i),'.pos');
end
% diff = w-u;
% namefile = sprintf('output/diff.pos');
% namesol = strcat('diff');
% writeField2D(dofm, mesh, diff, namefile, namesol);
% command = strcat(command, ' output/diff.pos');
% 
% [solP, sysP] = computeSolProjL2_2D_CG(mesh, dofm);
% L2norm = sqrt(solP'*Melim*solP);
% solP = solP./L2norm;
% namefile = sprintf('output/solP.pos');
% namesol = strcat('solP');
% writeField2D(dofm, mesh, solP, namefile, namesol);
% namefile = sprintf('output/rhsP.pos');
% namesol = strcat('rhsP');
% L2norm = sqrt(sysP.rhsP'*Melim*sysP.rhsP);
% sysP.rhsP = sysP.rhsP./L2norm;
% writeField2D(dofm, mesh, sysP.rhsP, namefile, namesol);
% command = strcat(command, ' output/solP.pos output/rhsP.pos');
% 
% 
% command = strcat(command, '&');
% system(command);


utMu = u'*Melim*u;
utu = uu'*uu;
vtMv = v'*Melim*v;
vtv = vv'*vv;
disp(['|| I - u^T*M*u || = ' num2str(norm(eye(nbEigVec)-utMu,2)) ]);
disp(['|| I - u^T*u || = ' num2str(norm(eye(nbEigVec)-utu,2)) ]);
disp(['|| I - v^T*M*v || = ' num2str(norm(eye(nbEigVec)-vtMv,2)) ]);
disp(['|| I - v^T*v || = ' num2str(norm(eye(nbEigVec)-vtv,2)) ]);


% utMu = u'*Melim*u;
% vtMv = v'*Melim*v;
% % utu = u'*u;
% % vtv = v'*v;
% % vptvp = vp'*v;
% u_v = u - v;
% disp(['|| I - ut*M*u || = ' num2str(norm(eye(nbEigVec)-utMu,2)) ]);
% disp(['|| I - vt*M*v || = ' num2str(norm(eye(nbEigVec)-vtMv,2)) ]);
% disp(['|| M || = ' num2str(norm(full(Melim),2)) ]);
% disp(['(u-v)t*M*(u-v) = ' num2str(sqrt(u_v'*Melim*u_v)) ]);


%% Compute eigs of A, M^-1*A, M0^-1*A and tabMin

mrange = 100; nrange = 100;
exactEigs = [];
for m = 1:mrange
    for n = 1:nrange
        exactEigs = [exactEigs; (m^2 + n^2)*pi^2];
    end
end
exactEigs = [exactEigs, exactEigs - k^2;];
exactEigs = sort(exactEigs,1);

% temp = sort(exactEigs(:,2),'ComparisonMethod','abs');


nbCompEigs = 50;
% vpA = eig(full(A));
% disp(['|  Spectrum of A computed']);
% vpMinvA = eig(full(A),full(M));
% disp(['|  Spectrum of M^{-1}A computed']);
% vpM0invA = eig(full(A),full(M0));
% disp(['|  Spectrum of M_{0}^{-1}A computed']);
vpA = eigs(A,nbCompEigs,'sm');
vpMinvA = eigs(A,M,nbCompEigs,'sm');
vpM0invA = eigs(A,M0,nbCompEigs,'sm');

vpA = sort(vpA);
vpMinvA = sort(vpMinvA);
vpM0invA = sort(vpM0invA);

spectrum = zeros(nbCompEigs,5);
spectrum(:, 1) = exactEigs(1:nbCompEigs,1);
spectrum(:, 2) = exactEigs(1:nbCompEigs,2);
spectrum(:,3) = vpA(1:nbCompEigs);
spectrum(:,4) = vpMinvA(1:nbCompEigs);
spectrum(:,5) = vpM0invA(1:nbCompEigs);

tabMin = zeros(nbCompEigs,5);
tabMin(:,1) = exactEigs(1:nbCompEigs,1);
tabMin(:,2) = exactEigs(1:nbCompEigs,2);
for i=1:nbCompEigs
    tabMin(i,3) = min(abs(vpA - tabMin(i,2)));
    tabMin(i,4) = min(abs(vpMinvA - tabMin(i,2)));
    tabMin(i,5) = min(abs(vpM0invA - tabMin(i,2)));
end

tabMinRel = tabMin;
tabMinRel(:,3:5) = tabMinRel(:,3:5)./tabMinRel(:,2);

figure;
hold on;
grid on;
plot(real(spectrum(1:nbCompEigs,2)), imag(spectrum(1:nbCompEigs,2)), 'o');
plot(real(spectrum(1:nbCompEigs,3)), imag(spectrum(1:nbCompEigs,3)), 's');
plot(real(spectrum(1:nbCompEigs,4)), imag(spectrum(1:nbCompEigs,4)), '^');
plot(real(spectrum(1:nbCompEigs,5)), imag(spectrum(1:nbCompEigs,5)), 'x');
leg = legend(['$\lambda_{m,n}^2 - k^2$'],['$A v = \lambda v$'], ['$M^{-1}A v = \lambda v$'], ['$(M_{0}^{-1}A) v = \lambda v$'], 'Location', 'best');
set(leg, 'Interpreter', 'latex', 'FontSize', 20);

%% Compute spectrum of deflated operators

tol = 1e-6;
k = 3.01*sqrt(2)*pi;
nbCompEigs = 10;
nbDefVec = 1;
eigsDef = zeros(nbCompEigs,3);

% Eigenvalies of M_0^-1*A
eigsDef(:,1) = eigs(Ared,Mred,nbCompEigs,'sm');

% Eigenvalues of Pdef*M_0^-1*A : construction of t<o different Pdef

% Deflation vectors
[u,nbDefVec,w] = computeProjEigVec_cavity_bis(mesh, dofm, nbDefVec,'closestEigvec',k);
u_int = R_int*u;
w_int = R_int*w;
u_int = u_int/sqrt((u_int)'*Mred*(u_int));
[v,~] = eigs(Ared,Mred,nbDefVec,'sm');
v = v/sqrt((v)'*Mred*(v));

% Deflation with proj
Z = u_int;
Zt = Z';
% MinvAZ = M\A*Z;
% E = Zt*MinvAZ;
% temp = E\Zt;
% EinvZtMinvA = temp/M*A;
% Q = Z*temp;
% Pdef = eye(size(A,1)) - M\A*Q;
% Qdef = eye(size(A,1)) - Q*M\A;
Pdef = eye(size(Ared,1)) - Mred\Ared* Z/(Zt*(Mred\Ared)*Z) * Zt;
PdefMinvA = (Mred\Ared) - (Mred\Ared) * Z/(Zt*(Mred\Ared)*Z) * Zt * (Mred\Ared);

% Deflation with gen eig vec
Z = v; Zt = Z';
Pdef_gen = eye(size(Ared,1)) - Mred\Ared* Z/(Zt*(Mred\Ared)*Z) * Zt;
PdefMinvA_gen = (Mred\Ared) - (Mred\Ared) * Z/(Zt*(Mred\Ared)*Z) * Zt * (Mred\Ared);

% Spectrum of Pdef*M_0^-1*A with proj
eigsDef(:,2) = eigs(PdefMinvA, nbCompEigs, 'sm');
% Spectrum of Pdef*M_0^-1*A with gen eigvec
eigsDef(:,3) = eigs(PdefMinvA_gen, nbCompEigs, 'sm');

figure;

% Ax=b
[uD, ~, ~, it, vecRes] = gmres(Ared, bred, size(Ared,1), tol, size(Ared,1));
disp(['    converged in ' num2str(size(vecRes,1)-1) ' iterations']);

vecRes = vecRes(:)./vecRes(1);

hold on;
semilogy(0:size(vecRes,1)-1, vecRes, 'DisplayName','Reduced : No def');

% M*Pdef*M^-1*Ax = M*Pdef*M^-1*b with proj vec
[uD, ~, ~, it, vecRes] = gmres(@(x) Mred*PdefMinvA*(x), Mred*(Pdef*(Mred\bred)), size(Ared,1), tol, size(Ared,1));
disp(['    converged in ' num2str(size(vecRes,1)-1) ' iterations']);

vecRes = vecRes(:)./vecRes(1);

hold on;
semilogy(0:size(vecRes,1)-1, vecRes, 'DisplayName','Reduced : Def with proj');

% M*Pdef*M^-1*Ax = M*Pdef*M^-1*b with gen vec
[uD, ~, ~, it, vecRes] = gmres(@(x) Mred*PdefMinvA_gen*(x), Mred*(Pdef_gen*(Mred\bred)), size(Ared,1), tol, size(Ared,1));
disp(['    converged in ' num2str(size(vecRes,1)-1) ' iterations']);

vecRes = vecRes(:)./vecRes(1);

hold on;
semilogy(0:size(vecRes,1)-1, vecRes, 'DisplayName','Reduced : Def with gen eigvec');



% Complete problem


% Ax=b
[uD, ~, ~, it, vecRes] = gmres(A, b, size(A,1), tol, size(A,1));
disp(['    converged in ' num2str(size(vecRes,1)-1) ' iterations']);
vecRes = vecRes(:)./vecRes(1);
hold on;
semilogy(0:size(vecRes,1)-1, vecRes, 'DisplayName','Original : No def', 'LineStyle','--')


% Pdef*Ax = Pdef*b with proj vec
[u,~] = computeProjEigVec_cavity_bis(mesh, dofm, nbDefVec,'closestEigvec',k);
u = u/sqrt((u)'*Melim*(u));
Z = u;
Zt = Z';
Pdef = eye(size(A,1)) - A* Z/(Zt*A*Z) * Zt;
PdefA = A - A * Z/(Zt*A*Z) * Zt * A;

[uD, ~, ~, it, vecRes] = gmres(@(x) PdefA*(x), Pdef*b, size(A,1), tol, size(A,1));
disp(['    converged in ' num2str(size(vecRes,1)-1) ' iterations']);
vecRes = vecRes(:)./vecRes(1);
hold on;
semilogy(0:size(vecRes,1)-1, vecRes, 'DisplayName','Original : Def with proj', 'LineStyle','--');

% Pdef*Ax = Pdef*b with gen vec
[v,~] = eigs(A,nbDefVec, 'sm');
v = v/sqrt((v)'*Melim*(v));
Z = v;
Zt = Z';
Pdef = eye(size(A,1)) - A* Z/(Zt*A*Z) * Zt;
PdefA = A - A * Z/(Zt*A*Z) * Zt * A;

[uD, ~, ~, it, vecRes] = gmres(@(x) PdefA*(x), Pdef*b, size(A,1), tol, size(A,1));
disp(['    converged in ' num2str(size(vecRes,1)-1) ' iterations']);
vecRes = vecRes(:)./vecRes(1);
hold on;
semilogy(0:size(vecRes,1)-1, vecRes, 'DisplayName','Original : Def with gen vec', 'LineStyle','--');


set(gca, 'YScale', 'log')
box on;
grid on;
ylim auto;
xlabel('Iteration number');
ylabel('Relative residual');
legend('Location', 'southwest');
title([num2str(nbDefVec) ' vectors deflated']);
ylim([tol 1]);
leg = legend(['$A_{II}x=b$'],['$M_{II}P_{def}M_{II}^{-1}A_{II} x = M_{II}PdefM_{II}^{-1}b$ with $u$'], ['$M_{II}P_{def}M_{II}^{-1}A_{II} x = M_{II}PdefM_{II}^{-1}b$ with $v$'], ['$Ax=b$'], ['$P_{def}Ax=P_{def}b$ with u'], ['$P_{def}Ax=P_{def}b$ with v'], 'Location', 'best');
set(leg, 'Interpreter', 'latex', 'FontSize', 12);

%% GMRES with right prec

[~,nbDefVec,w] = computeProjEigVec_cavity_bis(mesh, dofm, 1,'closestEigvec',k);
w_int = R_int*w;

figure;
% A*M^-1*x = b reduced

[uD, ~, ~, it, vecRes] = gmres(@(x) Ared*(Mred\(x)), bred, size(Ared,1), tol, size(Ared,1));
disp(['    converged in ' num2str(size(vecRes,1)-1) ' iterations']);

vecRes = vecRes(:)./vecRes(1);

hold on;
semilogy(0:size(vecRes,1)-1, vecRes);

% Pdef*A*M^-1*x = Pdef*b with M*proj vec
Z = w_int;
Zt = Z';
Pdef = eye(size(Ared,1)) - Ared/Mred* Z/(Zt*(Ared/Mred)*Z) * Zt;
PdefAMinv = (Ared/Mred) - (Ared/Mred) * Z/(Zt*(Ared/Mred)*Z) * Zt * Ared/Mred;

[uD, ~, ~, it, vecRes] = gmres(@(x) PdefAMinv*(x), Pdef*bred, size(Ared,1), tol, size(Ared,1));
disp(['    converged in ' num2str(size(vecRes,1)-1) ' iterations']);

vecRes = vecRes(:)./vecRes(1);

hold on;
semilogy(0:size(vecRes,1)-1, vecRes);


% A*M^-1*x = b 

[uD, ~, ~, it, vecRes] = gmres(@(x) A*(Melim\(x)), b, size(A,1), tol, size(A,1));
disp(['    converged in ' num2str(size(vecRes,1)-1) ' iterations']);

vecRes = vecRes(:)./vecRes(1);

hold on;
semilogy(0:size(vecRes,1)-1, vecRes, 'LineStyle','--');


set(gca, 'YScale', 'log')
box on;
grid on;
ylim auto;
xlabel('Iteration number');
ylabel('Relative residual');
legend('Location', 'southwest');
title([num2str(nbDefVec) ' vectors deflated']);
ylim([tol 1]);
leg = legend(['$A_{II}M_{II}^{-1}x=b$'],['$P_{def}A_{II}M_{II}^{-1} x = P_{def}b$ with $u$'], ['$AM_{elim}^{-1}x=b$'], 'Location', 'best');
set(leg, 'Interpreter', 'latex', 'FontSize', 12);



%% Batch freq
global k
k= 3.01*sqrt(2)*pi;
kmax = 13.8;
kmin = 12.7;
% rangeFreq = kmin:kmax;
rangeFreq = kmin:0.1:kmax;

nbit = zeros(length(rangeFreq),6);
nbit(:,1) = rangeFreq;

nbit(:,2) = run('red',rangeFreq,0,0,''); % Reducede system : Ax=b
nbit(:,3) = run('red',rangeFreq,1,3.01*sqrt(2)*pi,'proj'); % Reduced system : M*Pdef*M^-1*Ax = M*Pdef*M^-1*b with proj vec 
nbit(:,4) = run('red',rangeFreq,1,3.01*sqrt(2)*pi,'eigs'); % Reduced system : M*Pdef*M^-1*Ax = M*Pdef*M^-1*b with gen vec
nbit(:,5) = run('all',rangeFreq,0,0,''); % Original system : Ax=b
nbit(:,6) = run('all',rangeFreq,1,3.01*sqrt(2)*pi,'proj'); % Original system : PdefAx=Pdefb

figure;
hold on;
plot(nbit(:,1), nbit(:,2));
plot(nbit(:,1), nbit(:,3));
plot(nbit(:,1), nbit(:,4));
plot(nbit(:,1), nbit(:,5), 'LineStyle','--');
plot(nbit(:,1), nbit(:,6), 'LineStyle','--');

set(gca, 'YScale', 'log')
box on;
grid on;
ylim auto;
xlabel('Wavenumber');
ylabel('Number of iterations');
legend('Location', 'southwest');
title([num2str(nbDefVec) ' vectors deflated']);
leg = legend(['$A_{II}x=b$'],['$M_{II}P_{def}M_{II}^{-1}A_{II} x = M_{II}PdefM_{II}^{-1}b$ with $u$'], ['$M_{II}P_{def}M_{II}^{-1}A_{II} x = M_{II}PdefM_{II}^{-1}b$ with $v$'], ['$Ax=b$'], ['$P_{def}Ax=P_{def}b$ with u'], 'Location', 'best');
set(leg, 'Interpreter', 'latex', 'FontSize', 12);

function nbiter = run(domain,rangeFreq,nbEigVec,defFreq,type)
    global k;
    k = defFreq;
    mesh = setupBenchmark2D('cavity');
    mesh = buildConnectivity2D(mesh);
    dofm = buildDofManager2D_CG(mesh, 2);

    [~, sysA] = computeSolNum2D_CG_red(mesh, dofm, 0);

    switch domain
        case 'red'
            A = sysA.matAred;
            M = sysA.matMred;
            R_int = sysA.matR_int;
        case 'all'
            A = sysA.matA;
    end

    if nbEigVec > 0
        switch type
            case 'proj'
                [u,nbEigVec] = computeProjEigVec_cavity_bis(mesh, dofm, 1,'closestEigvec',defFreq);
                if (strcmp(domain,'red'))
                    u = R_int*u;
                end
            case 'eigs'
                [u,~] = eigs(A,M,1,'sm');
        end
        Z = u;
        Zt = Z';
    end

    nbiter = zeros(length(rangeFreq),1);
    i = 1;

    for k = rangeFreq
     
        [~, sysA] = computeSolNum2D_CG_red(mesh, dofm, 0);

        switch domain
            case 'red'
                A = sysA.matAred;
                M = sysA.matMred;
                b = sysA.rhsAred;
            case 'all'
                A = sysA.matA;
                b = sysA.rhsA;
        end

        if nbEigVec > 0
            switch domain
                case 'red'
                    op = M\A;
                case 'all'
                    op = A;
            end
            Pdef =  speye(size(A,1)) - op* Z/(Zt*op*Z) * Zt;
            PdefA = op - op * Z/(Zt*op*Z) * Zt * op;
            switch domain
                case 'red'
                    [~, ~, ~, it, ~] = gmres(@(x) M*PdefA*x,M*(Pdef*(M\b)),size(A,1),1e-6,size(A,1));
                case 'all'
                    [~, ~, ~, it, ~] = gmres(@(x) PdefA*x,Pdef*b,size(A,1),1e-6,size(A,1));
            end
        else
            [~, ~, ~, it, ~] = gmres(A,b,size(A,1),1e-6,size(A,1));
        end

    nbiter(i) = it(2) + (it(1)-1)*size(A,1);
    i = i + 1;


    end
end
