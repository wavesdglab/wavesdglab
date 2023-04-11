% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% GMRES with symmetric preconditioning

function [resVec, errorVec, i, flag, x] = solverGMRES(mesh, dofm, sys, tol, iMax, iOut, computeError)

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
    end
    %%%%%%%
    
    if (relRes <= tol)
        flag = 1;
        break;
    end
    i = i+1;
end
% Write the solution
%writeField2D(dofm, mesh, x, 'output/solNumGMRES.pos', "solNumGMRES");
%system('gmsh output/solNumGMRES.pos&');
end