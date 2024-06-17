% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% GMRES with symmetric preconditioning

function [resRedVec, resPhyVec, errorVec, i, flag, xPhy, x] = solverGMRESredu_DG_Marmousi(mesh, dofm, sys, sol, tol, restart, iMax, iOut, x0, it, computeError)

A = sys.matS;
b = sys.rhsS;
P = sys.matP;
Pinv = sys.matPinv;

x = x0;
count = 1;

H = zeros(restart+1,restart+1);
Q = zeros(size(A,2),restart+1);
sn = zeros(restart,1);
cs = zeros(restart,1);
beta = zeros(restart+1,1);

r = Pinv*(b-A*x);
beta(1) = sqrt(r'*P*r);
Q(:,1) = r/beta(1);

resRedVec = [];
resPhyVec = [];
errorVec  = [];

%%%%%%%

xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
rPhy = sys.rhsPhy - sys.matPhy*xPhy;
resRedIni = abs(beta(1));
resPhyIni = sqrt(rPhy'*rPhy);

if (~any(x))
    resRedVec(1) = 1;
    resPhyVec(1) = 1;
    errorVec(1) = computeError(mesh, dofm, xPhy, sol);
    fprintf('[%i] %g %g %g\n', 1, resRedVec(1), resPhyVec(1), errorVec(1));
end

%%%%%%%

flag = 0;
i = 1;
while(i <= restart)
    
%     if(mod(i,iOut) == 0)
%         [x,flag,relRes] = gmres(A,b,[],tol,i);
%         xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
%         rPhy = sys.rhsPhy - sys.matPhy*xPhy;
%         resPhyNew = rPhy'*rPhy;
%         resRedVec(i/iOut+1) = relRes;
%         resPhyVec(i/iOut+1) = sqrt(resPhyNew/resPhyIni);
%         errorVec(i/iOut+1) = computeNormError2D_DG(mesh, dofm, xPhy);
%         fprintf('[%i] %g %g\n', i, resRedVec(i/iOut+1), errorVec(i/iOut+1));
%     end
    
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
    relRes = abs(beta(i+1)) / resRedIni;

    %%%%%%%

    if(mod(i,restart) == 0)

        y = H(1:i,1:i) \ beta(1:i);
        x = Q(:,1:i) * y;

        if(mod((i+it*restart),iOut) == 0)
            computeError(mesh, dofm, xPhy, sol);
            xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
            rPhy = sys.rhsPhy - sys.matPhy*xPhy;
            resPhyNew = sqrt(rPhy'*rPhy);
            resRedVec(count) = relRes;
            resPhyVec(count) = resPhyNew/resPhyIni;
            errorVec(count) = computeError(mesh, dofm, xPhy, sol);
            fprintf('[%i] %g %g %g\n', i+it*restart, resRedVec(count), resPhyVec(count), errorVec(count));
            count = count + 1;
        end
    end
    %%%%%%%

    if (relRes <= tol)
        flag = 1;
        break;
    end
    i = i+1;
end

xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);

end