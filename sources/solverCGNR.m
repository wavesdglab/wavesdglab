% CGNR with symmetric preconditioning

function [resVec, errorVec, i, flag, x] = solverCGNR(mesh, dofm, sys, tol, iMax, iOut, computeError)

A = sys.matA;
b = sys.rhsA;
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

resVec = zeros(iMax/iOut+1,1);
errorVec = zeros(iMax/iOut+1,1);

%%%%%%%
resVec(1) = 1;
errorVec(1) = computeError(mesh, dofm, x);
fprintf('[%i] %g %g\n', 0, resVec(1), errorVec(1));
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
        resVec(i/iOut+1) = sqrt(rrnew/rrini);
        errorVec(i/iOut+1) = computeError(mesh, dofm, x);
        fprintf('[%i] %g %g\n', i, resVec(i/iOut+1), errorVec(i/iOut+1));
        %xRef = pcg(A'*A,A'*b,1e-10,i);
        %eRef = computeError(mesh, dofm, xRef);
        %fprintf('[%i] %g %g %g\n', i, resVec(i/iOut+1), errorVec(i/iOut+1), eRef);
    end
    %%%%%%%
    
    if (sqrt(rrnew/rrini) <= tol)
        flag = 1;
        break;
    end
    i = i+1;
end

end