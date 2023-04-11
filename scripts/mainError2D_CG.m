% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Authors: Axel Modave, Pierre Marchand, Timothée Raynaud

clear all;
%close all;

global k

% Setup benchmark and parameters
benchmark = 'cavity';
switch benchmark
    case 'open'
        k = 15*pi;
        h = 1/16;
        tol = 1e-10; maxit = 1000; itout = 50;
    case 'cavity'
        k = 0.1*sqrt(2)*pi;
        h = 1/16;
        tol = 1e-10; maxit = 2000; itout = 100;
    case 'waveguide'
        k = 6*pi;
        h = 1/8;
        tol = 1e-10; maxit = 4000; itout = 200;
end
degree = 1; % P1
PREC = 0; % for preconditioner

% Build mesh and DOF manager
mesh = benchmark2D(benchmark,h);
mesh = buildConnectivity2D(mesh);
dofm = buildDofManager2D_CG(mesh, degree); % espace fonctionnel discret

file = fopen('output/errorL2.dat','w');

fprintf(file,'%6s %12s %7s\n','k/(sqrt(2)*pi)','L2-Error (numSol)','Dlambda');

while(k<8.0*sqrt(2)*pi)

    Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1); % nb de points par longueur d'onde

    % -------------------------------------------------------------------------
    % Compute solution and error
    % -------------------------------------------------------------------------

    disp(['---------------------------------------------------------']);
    disp(['Method CG - Benchmark "' benchmark '"']);
    disp(['---------------------------------------------------------']);
    disp(['    k/(sqrt(2)*pi)     ' num2str(k/(sqrt(2)*pi))]);
    disp(['    h                   ' num2str(h)]);
    disp(['    degree              ' num2str(degree)]);
    disp(['    Dlambda             ' num2str(Dlambda)]);
    disp(['---------------------------------------------------------']);

    [solA, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);
    errorL2 = computeNormError2D_CG(mesh, dofm, solA);

    % solP = computeSolProjL2_2D_CG(mesh, dofm);
    % errorProjL2 = computeNormError2D_CG(mesh, dofm, solP);

    disp(['    L2-Error (numSol)   ' num2str(errorL2,'%1.2e')]);
    % disp(['    L2-Error (projSol)  ' num2str(errorProjL2,'%1.2e')]);
    disp(['---------------------------------------------------------']);

    fprintf(file,'%8.2f %17.2e %12.2f \n',k/(sqrt(2)*pi),errorL2,Dlambda);
    
    k = k + 0.1*sqrt(2)*pi;

end
fclose(file);
