% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% GMRES with symmetric preconditioning

function [resVec, errorVec, i, flag, x] = solverGMRES(mesh, dofm, sys, tol, iMax, iOut, computeError)


% [solA, ~] = computeSolNum2D_CG(mesh, dofm, 1);
% errorL2 = computeNormError2D_CG(mesh, dofm, solA);

% solP = computeSolProjL2_2D_CG(mesh, dofm);
% errorProjL2 = computeNormError2D_CG(mesh, dofm, solP);


A = sys.matA;
b = sys.rhsA;
P = sys.matP;
Pinv = sys.matPinv;

x = zeros(size(A,2),1);
H = zeros(iMax+1,iMax+1);
Q = zeros(size(A,2),iMax+1);
sn = zeros(iMax,1);
cs = zeros(iMax,1);
beta = zeros(iMax+1,1);

r = Pinv*(b-A*x);
beta(1) = sqrt(r'*P*r);
Q(:,1) = r/beta(1);

resVec = zeros(iMax/iOut+1,1);
errorVec = zeros(iMax/iOut+1,1);

%%%%%%%
resVec(1) = 1;
errorVec(1) = computeError(mesh, dofm, x);
% fprintf('[%i] %g %g\n', 0, resVec(1), errorVec(1));
%%%%%%%
%%%%%%%
eigenvalArray = zeros(iMax, size(A,2));

% figure
% set(0,'DefaultFigureWindowStyle','docked')

% hold on

% subplot(2,1,1);
% mat = Pinv * A;
% [eigenvec, eigenval] = eigs(mat,size(mat,1));
% eigenvalA = diag(eigenval);
% scatter(real(eigenvalA),imag(eigenvalA), 'r', 'Marker', 'x');
% xlim([-200 200]);
% ylim([-1 1]);
% % title(['Harmonic Ritz values at iteration 0'], 'interpreter', 'latex', 'fontsize', 20)
% grid on; box on;

% subplot(2,1,2);
% plot([0 iMax],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (direct)','linewidth', 1);
% plot([0 iMax],[errorProjL2 errorProjL2],'k:','DisplayName','Relative L2-error (projection)','linewidth', 1);
% semilogy(0,resVec(1),'b-o','DisplayName','Relative residual','linewidth', 1,'markersize', 5);
% semilogy(0,errorVec(1),'k-o','DisplayName','Relative L2-error (iterative)','linewidth', 1,'markersize', 5);
% set(gca, 'YScale', 'log')
% box on
% grid on
% xlim([0 155]);
% ylim([10e-12 10e2]);
% title(['CG - Cavity - GMRES - k=2*\sqrt(2)*\pi - h=1/16 - degree=1'], 'interpreter', 'latex', 'fontsize', 20)
% xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
% ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
% legend('Location', 'southwest', 'fontsize', 15)
% grid on; box on;


%%%%%%%

flag = 0;
i = 1;
while(i <= iMax)
    
    % Arnoldi iteration – Add one vector to basis Q and orthogonalize it
    Q(:,i+1) = Pinv*A*Q(:,i);
    for j = 1:i
        H(j,i) = Q(:,j)' * P * Q(:,i+1);
        Q(:,i+1) = Q(:,i+1) - H(j,i) * Q(:,j);
    end
    H(i+1,i) = sqrt(Q(:,i+1)' * P * Q(:,i+1));
    Q(:,i+1) = Q(:,i+1) / H(i+1,i);

    % Apply the previous Givens matrix to ith column
    for j = 1:i-1
        matGivens = [ cs(j)' sn(j)' ; -sn(j) cs(j) ];
        H(j:j+1,i) = matGivens * H(j:j+1,i);
    end
    
    % Compute the new Givens matrix
    tmp = sqrt(abs(H(i,i))^2 + H(i+1,i)^2);
    cs(i) = H(i,i)/tmp;    % complex
    sn(i) = H(i+1,i)/tmp;  % real
    matGivens = [ cs(i)' sn(i)' ; -sn(i) cs(i) ];
    
    % Apply the new Givens matrix to ith column of H and residual vector
    H(i:i+1,i)  = matGivens * H(i:i+1,i);
    beta(i:i+1) = matGivens * beta(i:i+1);
    
    % Update the residual vector
    relRes = abs(beta(i+1)) / norm(beta(1));
    
    %%%%%%%
    if(mod(i,iOut) == 0)
        y = H(1:i,1:i) \ beta(1:i);
        x = Q(:,1:i) * y;
 
        resVec(i/iOut+1) = relRes;
        errorVec(i/iOut+1) = computeError(mesh, dofm, x);
        fprintf('[%i] %g %g\n', i, resVec(i/iOut+1), errorVec(i/iOut+1));
        %xRef = gmres(A,b,[],1e-10,i);
        %eRef = computeError(mesh, dofm, xRef);
        %fprintf('[%i] %g %g %g\n', i, resVec(i/iOut+1), errorVec(i/iOut+1), eRef);

        matQ = eye(i,i);
        for j = 1:i-1
            tempQ = eye(i,i);
            tempQ(j,j) = cs(j)';
            tempQ(j+1,j+1) = cs(j);
            tempQ(j,j+1) = sn(j)';
            tempQ(j+1,j) = -sn(j)';
            matQ = tempQ*matQ;
        end

        H_sub = zeros(i, i);


        for j = 1:i
            H_sub(j, 1:i) = H(j, 1:i);
        end

        [V,D] = eig(H_sub'*H_sub,H_sub'*matQ);
        eigenval = diag(D);

        eigenvalArray(i, 1:length(eigenval)) = eigenval';

        % clf;
        % set(0,'DefaultFigureWindowStyle','docked')
        % subplot(2,1,1);
        % grid on; box on;
        

        % hold on
        % scatter(real(eigenvalA),imag(eigenvalA), 'r', 'Marker', 'x');
        % scatter(real(eigenval),imag(eigenval), 'b');
        % xlim([-200 200]);
        % ylim([-1 1]);
        
        
        % title(['Harmonic Ritz values at iteration ' num2str(i)], 'interpreter', 'latex', 'fontsize', 20)
        % drawnow;
        
        % nbOut = floor(i/iOut) +1;

        % subplot(2,1,2);
        % hold on
        % plot([0 iMax],[errorL2 errorL2],'k--','DisplayName','Relative L2-error (direct)','linewidth', 1);
        % plot([0 iMax],[errorProjL2 errorProjL2],'k:','DisplayName','Relative L2-error (projection)','linewidth', 1);
        % semilogy(0:iOut:i,resVec(1:nbOut),'b-o','DisplayName','Relative residual','linewidth', 1,'markersize', 5);
        % semilogy(0:iOut:i,errorVec(1:nbOut),'k-o','DisplayName','Relative L2-error (iterative)','linewidth', 1,'markersize', 5);
        % set(gca, 'YScale', 'log')
        % box on
        % grid on
        % xlim([0 155]);
        % ylim([10e-12 10e2]);
        % title(['CG - Cavity - GMRES - k= ' num2str(2*sqrt(2)*pi) ' - h=' num2str(1/16) ' - degree=' num2str(1)], 'interpreter', 'latex', 'fontsize', 20)
        % xlabel('Iteration', 'interpreter', 'Latex', 'fontsize', 15)
        % ylabel('Values', 'interpreter', 'Latex', 'fontsize', 15)
        % legend('Location', 'southwest', 'fontsize', 15)
        % grid on; box on;
        % drawnow;

        % pause(0.5);
    end
    %%%%%%%
    
    if (relRes <= tol)
        flag = 1;
        break;
    end
    i = i+1;
end
csvwrite('output/hrv.csv', eigenvalArray);
csvwrite('output/res_rel.csv', resVec);
% Write the solution
%writeField2D(dofm, mesh, x, 'output/solNumGMRES.pos', "solNumGMRES");
%system('gmsh output/solNumGMRES.pos&');
end