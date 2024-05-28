% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Preconditioned CGNR

function [resRedVec, resPhyVec, errorVec, i, flag, xPhy] = solverCGNRredu_DG(mesh, dofm, sys, tol, iMax, iOut, computeError)

A = sys.matS;
b = sys.rhsS;
P = sys.matP;
Pinv = sys.matPinv;

% % Right-preconditioning
% x = zeros(size(A,2),1);
% r = b - A*x;
% z = Pinv*(A'*r);
% p = z;
% rrini = r'*r;
% zzini = z'*z;
% zzold = zzini;

% Symmetric-preconditioning
x = zeros(size(A,2),1);
r = b - A*x;
z = Pinv*(A'*(Pinv*r));
p = z;
rrini = r'*r;
zzini = z'*P*z;
zzold = zzini;

resRedVec = zeros(iMax/iOut+1,1);
resPhyVec = zeros(iMax/iOut+1,1);
errorVec  = zeros(iMax/iOut+1,1);

%%%%%%%
xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
rPhy = sys.rhsPhy - sys.matPhy*xPhy;
resPhyIni = sqrt(rPhy'*rPhy);
resRedVec(1) = 1;
resPhyVec(1) = 1;
errorVec(1) = computeError(mesh, dofm, xPhy);
fprintf('[%i] %g %g %g\n', 1, resRedVec(1), resPhyVec(1), errorVec(1));
%%%%%%%

flag = 0;
i = 1;
while(i <= iMax)
    
%     % Right-preconditioning
%     q = A*(Pinv*p);
%     alpha = zzold/(q'*q);
%     x = x + alpha*(Pinv*p);
%     r = r - alpha*q;
%     z = Pinv*(A'*r);
%     rrnew = r'*r;
%     zznew = z'*z;
%     beta = zznew/zzold;
%     p = z + beta*p;
%     zzold = zznew;
    
    % Symmetric-preconditioning
    q = Pinv*(A*p);
    alpha = zzold/(q'*P*q);
    x = x + alpha*p;
    r = r - alpha*A*p;
    z = Pinv*(A'*(Pinv*r));
    rrnew = r'*r;
    zznew = z'*P*z;
    beta = zznew/zzold;
    p = z + beta*p;
    zzold = zznew;
    
    %%%%%%%
    if(mod(i,iOut) == 0)
        xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
        rPhy = sys.rhsPhy - sys.matPhy*xPhy;
        resPhyNew = sqrt(rPhy'*rPhy);
        resRedVec(i/iOut+1) = sqrt(rrnew/rrini);
        resPhyVec(i/iOut+1) = resPhyNew/resPhyIni;
        errorVec(i/iOut+1) = computeError(mesh, dofm, xPhy);
        fprintf('[%i] %g %g %g\n', i, resRedVec(i/iOut+1), resPhyVec(i/iOut+1), errorVec(i/iOut+1));
    end
    %%%%%%%
    
    if (sqrt(rrnew/rrini) <= tol)
        flag = 1;
        break;
    end
    i = i+1;
end

xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);

end