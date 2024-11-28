% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% GMRES with symmetric preconditioning

function [resRedVec, resPhyVec, errorVec, i, flag, xPhy] = solverGMRESredu_DG(mesh, dofm, sys, tol, iMax, iOut, computeError, xRef)

A = sys.matS;
b = sys.rhsS;
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

resRedVec = zeros(iMax/iOut+1,1);
resPhyVec = zeros(iMax/iOut+1,1);
errorVec  = zeros(iMax/iOut+1,1);

%%%%%%%
xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
rPhy = sys.rhsPhy - sys.matPhy*xPhy;
resRedIni = abs(beta(1));
resPhyIni = sqrt(rPhy'*rPhy);
resRedVec(1) = 1;
resPhyVec(1) = 1;
errorVec(1) = computeError(mesh, dofm, xPhy, xRef);
fprintf('[%i] %g %g %g\n', 1, resRedVec(1), resPhyVec(1), errorVec(1));
%%%%%%%

flag = 0;
i = 1;
while(i <= iMax)
    
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
    if(mod(i,iOut) == 0)
        y = H(1:i,1:i) \ beta(1:i);
        x = Q(:,1:i) * y;
        
        xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
        rPhy = sys.rhsPhy - sys.matPhy*xPhy;
        resPhyNew = sqrt(rPhy'*rPhy);
        resRedVec(i/iOut+1) = relRes;
        resPhyVec(i/iOut+1) = resPhyNew/resPhyIni;
        errorVec(i/iOut+1) = computeError(mesh, dofm, xPhy, xRef);
        fprintf('[%i] %g %g %g\n', i, resRedVec(i/iOut+1), resPhyVec(i/iOut+1), errorVec(i/iOut+1));
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