% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% CGNR with symmetric preconditioning

function [resRedVec, resPhyVec, errorVec, i, flag, xPhy] = solverCGNRredu_CG(mesh, dofm, sys, tol, iMax, iOut, computeError)

A = sys.matS;
b = sys.rhsS;
Pinv = sys.matPinv;

x = zeros(size(A,2),1);
r = b-A*x;
s = Pinv*r;  % s=r for left-preconditioning
y = A'*s;
z = Pinv*y;
p = z;
rrini = r'*r;
zzini = y'*z;
zzold = zzini;

resRedVec = zeros(iMax/iOut+1,1);
resPhyVec = zeros(iMax/iOut+1,1);
errorVec  = zeros(iMax/iOut+1,1);

%%%%%%%
xPhy = [ x ; sys.matIIinv*(sys.rhsI-sys.matIG*x) ];
rPhy = sys.rhsA - sys.matA*xPhy;
resPhyIni = rPhy'*rPhy;
resRedVec(1) = 1;
resPhyVec(1) = 1;
errorVec(1) = computeError(mesh, dofm, xPhy);
fprintf('[%i] %g %g\n', 0, resRedVec(1), errorVec(1));
%%%%%%%

flag = 0;
i = 1;
while(i <= iMax)
    
    v = A*p;
    w = Pinv*v;  % w=v for left-preconditioning
    alpha = zzold/(v'*w);
    x = x + alpha*p;
    r = r - alpha*v;
    s = Pinv*r;  % s=r for left-preconditioning
    y = A'*s;
    z = Pinv*y;
    rrnew = r'*r;
    zznew = y'*z;
    p = z + (zznew/zzold)*p;
    zzold = zznew;
    
    %%%%%%%
    if(mod(i,iOut) == 0)
        xPhy = [ x ; sys.matIIinv*(sys.rhsI-sys.matIG*x) ];
        rPhy = sys.rhsA - sys.matA*xPhy;
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

xPhy = [ x ; sys.matIIinv*(sys.rhsI-sys.matIG*x) ];

end