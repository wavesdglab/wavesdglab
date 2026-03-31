%% Setup
%close all;
clear;
clear global;

global k h Options

Options.Basis = 'Jacobi'; % Jacobi, Lagrange
Options.Error = 'L2'; % L2, H1


%%%%%%%%%%%%%%%%%%%%%%%% Open Cavity

global LdomX LdomY LpmlX LpmlY WRITE_FIELD_ABSOLUTE
LdomX = 0.95; LdomY = 0.5; LpmlX = 0.2; LpmlY = 0.2;
WRITE_FIELD_ABSOLUTE = 0;

k = 23.591;
h = 1/10;
tol = 1e-4;
degree = 3;

mesh = setupBenchmark2D('scattering_openCavity_NEU'); % cavity or scattering_openCavity_NEU
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree);
Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

[~, sysA] = computeSolNum2D_CG_red(mesh, dofm, 'none');
A = sysA.matA;
P = sysA.matP;
M = sysA.matM;
b = sysA.rhsA;
I = eye(size(A,2));

M0 = sysA.matMevpb;
Melim = sysA.matMelim;
Mred = sysA.matMred;
Ared = sysA.matAred;
bred = sysA.rhsAred;

R_int = sysA.matR_int;
R = computeRestrictionCav(mesh,dofm);

%% M-orthogonality

nbEigVec = 3;
[u,nbEigVec] = computeProjEigVec_openCavity_NEU_bis(mesh, dofm, nbEigVec, k);
[v,~] = eigs(A,M0,nbEigVec,'sm');
% [v,~] = eigs(A,nbEigVec,'sm');
% u = R*u;
% v = R*v;
uu = u;
vv = v;
command = strcat('gmsh output/mesh.msh');
for i = 1:nbEigVec
    L2norm = sqrt(abs(u(:,i)'*Melim*u(:,i)));
    u(:,i) = -u(:,i)./L2norm;
    namefile = sprintf('output/projvec_%g.pos',i);
    namesol = strcat('projvec_',num2str(i));
    writeField2D(dofm, mesh, u(:,i), namefile, namesol);

    L2norm = sqrt(abs(v(:,i)'*Melim*v(:,i)));
    v(:,i) = v(:,i)./L2norm;
    namefile = sprintf('output/eigvec_%g.pos',i);
    namesol = strcat('eigvec_',num2str(i));
    writeField2D(dofm, mesh, v(:,i), namefile, namesol);
    command = strcat(command, ' output/projvec_',num2str(i),'.pos output/eigvec_',num2str(i),'.pos');

    L2norm = sqrt(abs(uu(:,i)'*uu(:,i)));
    uu(:,i) = uu(:,i)./L2norm;
    L2norm = sqrt(abs(vv(:,i)'*vv(:,i)));
    vv(:,i) = vv(:,i)./L2norm;
end
diff = v-u;
% namefile = sprintf('output/diff.pos');
% namesol = strcat('diff');
% writeField2D(dofm, mesh, diff, namefile, namesol);
% command = strcat(command, ' output/diff.pos');
% 
command = strcat(command, '&');
system(command);



utMu = u'*Melim*u;
vtMv = v'*Melim*v;
utu = uu'*uu;
vtv = vv'*vv;
disp(['|| I - u^T*M*u || = ' num2str(norm(eye(nbEigVec)-utMu,2)) ]);
disp(['|| I - u^T*u || = ' num2str(norm(eye(nbEigVec)-utu,2)) ]);
disp(['|| I - v^T*M*v || = ' num2str(norm(eye(nbEigVec)-vtMv,2)) ]);
disp(['|| I - v^T*v || = ' num2str(norm(eye(nbEigVec)-vtv,2)) ]);
disp(['(u-v)t*M*(u-v) = ' num2str(sqrt(diff(:,1)'*Melim*diff(:,1))) ]);


%% Compute spectrum of deflated operators

nbCompEigs = 20;
nbDefVec = 1;
eigsDef = zeros(nbCompEigs,3);

% Eigenvalies of M_0^-1*A
eigsDef(:,1) = eigs(Ared,Mred,nbCompEigs,'sm');

% Eigenvalues of Pdef*M_0^-1*A : construction of t<o different Pdef

% Deflation vectors
[u,nbDefVec] = computeProjEigVec_openCavity_NEU_bis(mesh, dofm, nbDefVec, k);
u_int = R_int*u;
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
[u,~] = computeProjEigVec_openCavity_NEU_bis(mesh, dofm, nbDefVec, k);
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
