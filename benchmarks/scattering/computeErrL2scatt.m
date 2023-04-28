% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Authors: Axel Modave, Pierre Marchand, Timothée Raynaud

close all;
clear all;

global k R_disk L L_PML;

k = 25;
L = 1.1;
R_disk = 1;

degree = 1;
PREC = 0;
tol = 1e-10; maxit = 1000; itout = 50;



benchmark = 'scatteringPML';

L_PML_tab = [0.1, 0.2, 0.4];

h_tab = [0.1, 0.05, 0.025, 0.0125];

file = fopen('output/errL2.csv','w');

% pour chaque benchmark, pour chaque h, on calcule l'erreur en norme L2

for i = 1:length(L_PML_tab)
    L_PML = L_PML_tab(i);
    for j = 1:length(h_tab)
        h = h_tab(j);
        mesh = benchmark2D(benchmark,h);
        mesh = buildConnectivity2D(mesh);
        dofm = buildDofManager2D_CG(mesh, degree);

        Dlambda = 2*pi/k * (sqrt(dofm.numDofTRI) - 1);

        disp(['---------------------------------------------------------']);
        disp(['Method CG - Benchmark "' benchmark '"']);
        disp(['---------------------------------------------------------']);
        disp(['    k                   ' num2str(k)]);
        disp(['    h                   ' num2str(h)]);
        disp(['    degree              ' num2str(degree)]);
        disp(['    Dlambda             ' num2str(Dlambda)]);
        disp(['---------------------------------------------------------']);


        [solA, sysA] = computeSolNum2D_CG(mesh, dofm, PREC);
        errornumL2 = computeNormError2D_CG(mesh, dofm, solA);
        disp(['    L2-Error (numSol)   ' num2str(errornumL2,'%1.2e')]);

        solP = computeSolProjL2_2D_CG(mesh, dofm);
        errorL2 = computeNormError2D_CG(mesh, dofm, solP);
        disp(['    L2-Error (refSol)   ' num2str(errorL2,'%1.2e')]);


        if(L_PML == 0.4)
            solP = computeSolProjL2_2D_CG(mesh, dofm);
            errorL2 = computeNormError2D_CG(mesh, dofm, solP);
            disp(['    L2-Error (refSol)   ' num2str(errorL2,'%1.2e')]);

            fprintf(file, '%s, %f, %f, %f\n', benchmark, h, errornumL2, errorL2);
        else
            fprintf(file, '%s, %f, %f\n', benchmark, h, errornumL2);
        end
    end
end
fclose(file);