% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% CGNE with symmetric preconditioning

function [resRedVec, resPhyVec, errorVec, i, flag, xPhy] = solverCGNEredu_DG(mesh, dofm, sys, tol, iMax, iOut, computeError)

A = sys.matS;
b = sys.rhsS;
Pinv = sys.matPinv;

x = zeros(size(A,2),1);
r = b - A*x;
s = Pinv*r;
p = A'*s;
q = Pinv*p;
rr = r'*s;
rrini = rr;

resRedVec = zeros(iMax/iOut+1,1);
resPhyVec = zeros(iMax/iOut+1,1);
errorVec  = zeros(iMax/iOut+1,1);

%%%%%%%
xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
rPhy = sys.rhsPhy - sys.matPhy*xPhy;
resPhyIni = rPhy'*rPhy;
resRedVec(1) = 1;
resPhyVec(1) = 1;
errorVec(1) = computeError(mesh, dofm, xPhy);
fprintf('[%i] %g %g\n', 0, resRedVec(1), errorVec(1));
%%%%%%%

flag = 0;
i = 1;
while(i <= iMax)
    
    pp = p'*q;
    alpha = rr/pp;
    x = x + alpha*q;
    r = r - alpha*A*q;
    s = Pinv*r;
    rrnew = r'*s;
    beta = rrnew/rr;
    rr = rrnew;
    p = A'*s + beta*p;
    q = Pinv*p;
    
    %%%%%%%
    if(mod(i,iOut) == 0)
        xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
        rPhy = sys.rhsPhy - sys.matPhy*xPhy;
        resPhyNew = rPhy'*rPhy;
        resRedVec(i/iOut+1) = sqrt(rrnew/rrini);
        resPhyVec(i/iOut+1) = sqrt(resPhyNew/resPhyIni);
        errorVec(i/iOut+1) = computeError(mesh, dofm, xPhy);
        fprintf('[%i] %g %g\n', i, resRedVec(i/iOut+1), errorVec(i/iOut+1));
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