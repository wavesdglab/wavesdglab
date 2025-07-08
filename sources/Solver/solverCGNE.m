% CGNE with symmetric preconditioning

function [resVec, errorVec, i, flag, x] = solverCGNE(mesh, dofm, sys, tol, iMax, iOut, computeError)

A = sys.matA;
b = sys.rhsA;
P = sys.matP;

x = zeros(size(A,2),1);
r = b - A*x;
s = P\r;
p = A'*s;
q = P\p;
rr = real(r'*s);
rrini = rr;

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
    
    pp = real(p'*q);
    alpha = rr/pp;
    x = x + alpha*q;
    r = r - alpha*A*q;
    s = P\r;
    rrnew = real(r'*s);
    beta = rrnew/rr;
    rr = rrnew;
    p = A'*s + beta*p;
    q = P\p;
    
    %%%%%%%
    if(mod(i,iOut) == 0)
        resVec(i/iOut+1) = sqrt(rrnew/rrini);
        errorVec(i/iOut+1) = computeError(mesh, dofm, x);
        fprintf('[%i] %g %g\n', i, resVec(i/iOut+1), errorVec(i/iOut+1));
    end
    %%%%%%%
    
    if (sqrt(rrnew/rrini) <= tol)
        flag = 1;
        break;
    end
    i = i+1;
end

end