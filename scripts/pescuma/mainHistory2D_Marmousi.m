% close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MARMOUSI
iMax = 1000; iOut = 50; restart = 10;
benchmark = 'geophysics_marmousi';

degree = 3; freq = 5; tol = 1e-100;
run(benchmark,degree,tol,iMax,iOut,restart,freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(benchmark,degree,tol,iMax,iOut,restart,freq)
global omega nLambda
global rhoArray cArray

omega = 2*pi*freq;
nLambda = 10/(degree+1);

% Build mesh and DOF manager
mesh = setupBenchmark2D(benchmark);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_DG(mesh, degree);

% Print coefficients
writeCoef2D(mesh, cArray, 'output/velocity.pos', "Velocity [m/s]");
writeCoef2D(mesh, rhoArray, 'output/density.pos', "Density");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------------------------------------------');
disp(['Method CHDG (' benchmark ')']);
disp('---------------------------------------------------------');
disp(['    degree              ' num2str(degree)]);
disp('---------------------------------------------------------');

PREC = 1;

[solA, sysA] = computeSolNum2D_CHDG_heterogeneous_2(mesh, dofm, PREC);
[solB, sysB] = computeSolNum2D_HDG_heterogeneous_2(mesh, dofm, PREC); 

disp('---------------------------------------------------------');

% writeField2D(dofm, mesh, solA, 'output/solNum.pos', "solNum");
% system('gmsh output/solNum.pos&');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--- Solver CGNR');
[resRedVec1A, resPhyVec1A, error1A] = solverCGNRredu_DG_Marmousi(mesh, dofm, sysA, solA, tol, iMax, iOut, @computeNormError2D_DG_Marmousi);
[resRedVec1B, resPhyVec1B, error1B] = solverCGNRredu_DG_Marmousi(mesh, dofm, sysB, solB, tol, iMax, iOut, @computeNormError2D_DG_Marmousi);

iterVec = (0:iOut:iMax)';

rezu1A = ["iter" "resRed" "resPhy" "error"];
rezu2A = [iterVec resRedVec1A, resPhyVec1A, error1A];
name = sprintf('output/historyCGNR_CHDG_0thorder_%s_p%i_omega%g.csv', benchmark, degree, omega);
writematrix([rezu1A ; rezu2A], name, 'Delimiter', 'semi');

rezu1B = ["iter" "resRed" "resPhy" "error"];
rezu2B = [iterVec resRedVec1B, resPhyVec1B, error1B];
name = sprintf('output/historyCGNR_HDG_0thorder_%s_p%i_omega%g.csv', benchmark, degree, omega);
writematrix([rezu1B ; rezu2B], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--- Solver GMRES');

% First step of GMRES with restart, using a null initial guess x0
A = sysA.matS;
x0 = zeros(size(A,2),1);
it = 0;
[resRedVec3Anew, resPhyVec3Anew, error3Anew, ~, ~, ~, xS] = solverGMRESredu_DG_Marmousi(mesh, dofm, sysA, solA, tol, restart, iMax, iOut, x0, it, @computeNormError2D_DG_Marmousi);
resRedVec3Anew = resRedVec3Anew';
resPhyVec3Anew = resPhyVec3Anew';
error3Anew = error3Anew';

% Update initial guess for second step of GMRES with restart
x0 = xS;

T = iMax / restart;

for it = 1:T-1

    [resRedVec3A, resPhyVec3A, error3A, ~, ~, ~, xS] = solverGMRESredu_DG_Marmousi(mesh, dofm, sysA, solA, tol, restart, iMax, iOut, x0, it, @computeNormError2D_DG_Marmousi);

    resRedVec3Anew = [resRedVec3Anew; resRedVec3A'];
    resPhyVec3Anew = [resPhyVec3Anew; resPhyVec3A'];
    error3Anew = [error3Anew; error3A'];

    % Update initial guess for restart #(it+1)
    x0 = xS;

end

resRedVec3A = resRedVec3Anew;
resPhyVec3A = resPhyVec3Anew;
error3A = error3Anew;

iterVec = (0:iOut:iMax)';

rezu1A = ["iter" "resRed" "resPhy" "error"];
rezu2A = [iterVec resRedVec3A, resPhyVec3A, error3A];
name = sprintf('output/historyGMRES_CHDG_0thorder_%s_p%i_omega%g.csv', benchmark, degree, omega);
writematrix([rezu1A ; rezu2A], name, 'Delimiter', 'semi');



% First step of GMRES with restart, using a null initial guess x0
A = sysB.matS;
x0 = zeros(size(A,2),1);
it = 0;
[resRedVec3Bnew, resPhyVec3Bnew, error3Bnew, ~, ~, ~, xS] = solverGMRESredu_DG_Marmousi(mesh, dofm, sysB, solB, tol, restart, iMax, iOut, x0, it, @computeNormError2D_DG_Marmousi);
resRedVec3Bnew = resRedVec3Bnew';
resPhyVec3Bnew = resPhyVec3Bnew';
error3Bnew = error3Bnew';

% Update initial guess for second step of GMRES with restart
x0 = xS;

T = iMax / restart;

for it = 1:T-1

    [resRedVec3B, resPhyVec3B, error3B, ~, ~, ~, xS] = solverGMRESredu_DG_Marmousi(mesh, dofm, sysB, solB, tol, restart, iMax, iOut, x0, it, @computeNormError2D_DG_Marmousi);

    resRedVec3Bnew = [resRedVec3Bnew; resRedVec3B'];
    resPhyVec3Bnew = [resPhyVec3Bnew; resPhyVec3B'];
    error3Bnew = [error3Bnew; error3B'];

    % Update initial guess for restart #(it+1)
    x0 = xS;

end

resRedVec3B = resRedVec3Bnew;
resPhyVec3B = resPhyVec3Bnew;
error3B = error3Bnew;

iterVec = (0:iOut:iMax)';

rezu1B = ["iter" "resRed" "resPhy" "error"];
rezu2B = [iterVec resRedVec3B, resPhyVec3B, error3B];
name = sprintf('output/historyGMRES_HDG_0thorder_%s_p%i_omega%g.csv', benchmark, degree, omega);
writematrix([rezu1B ; rezu2B], name, 'Delimiter', 'semi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold off
semilogy(iterVec,resPhyVec1A,'-ob','MarkerFaceColor','b','DisplayName','CHDG (0th-order) - CGNR');
hold on
semilogy(iterVec,resPhyVec3A,'-ob','MarkerFaceColor','w','DisplayName','CHDG (0th-order) - GMRES');
semilogy(iterVec,resPhyVec1B,'-or','MarkerFaceColor','r','DisplayName','HDG - CGNR');
semilogy(iterVec,resPhyVec3B,'-or','MarkerFaceColor','w','DisplayName','HDG - GMRES');
box on;
grid on;
legend('Location','southwest');
xlabel('Iteration');
ylabel('Norm physical residual');
axis([0 iMax 1e-8 1]);

figure;
hold off
semilogy(iterVec,error1A,'-ob','MarkerFaceColor','b','DisplayName','CHDG (0th-order) - CGNR');
hold on
semilogy(iterVec,error3A,'-ob','MarkerFaceColor','w','DisplayName','CHDG (0th-order) - GMRES');
semilogy(iterVec,error1B,'-or','MarkerFaceColor','r','DisplayName','HDG - CGNR');
semilogy(iterVec,error3B,'-or','MarkerFaceColor','w','DisplayName','HDG - GMRES');
box on;
grid on;
legend('Location','southwest');
xlabel('Iteration');
ylabel('Relative error');
axis([0 iMax 1e-8 1]);

figure;
hold off
semilogy(iterVec,resRedVec1A,'-ob','MarkerFaceColor','b','DisplayName','CHDG (0th-order) - CGNR');
hold on
semilogy(iterVec,resRedVec3A,'-ob','MarkerFaceColor','w','DisplayName','CHDG (0th-order) - GMRES');
semilogy(iterVec,resRedVec1B,'-or','MarkerFaceColor','r','DisplayName','HDG - CGNR');
semilogy(iterVec,resRedVec3B,'-or','MarkerFaceColor','w','DisplayName','HDG - GMRES');
box on;
grid on;
legend('Location','southwest');
xlabel('Iteration');
ylabel('Norm reduced residual');
axis([0 iMax 1e-8 1]);

end