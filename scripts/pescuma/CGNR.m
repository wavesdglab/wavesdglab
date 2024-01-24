% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% CGNR with symmetric preconditioning

function [resRedVec] = CGNR(sys, tol, iMax, iOut)

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

%%%%%%%
resRedVec(1) = 1;
fprintf('[%i] %g\n', 0, resRedVec(1));
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
        resRedVec(i/iOut+1) = sqrt(rrnew/rrini);
        fprintf('[%i] %g\n', i, resRedVec(i/iOut+1));
    end
    %%%%%%%
    
    if (sqrt(rrnew/rrini) <= tol)
        flag = 1;
        break;
    end

    i = i+1;
end

end