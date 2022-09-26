function [resVec, errorVec, i, flag] = solverCGN_DG(mesh, dofm, sys, tol, iMax, iOut)

A = sys.matA;
b = sys.rhsA;
x = zeros(size(A,2),1);
r = b - A*x;

z = A'*r;
p = z;
rrini = r'*r;
zzini = z'*z;
zzold = zzini;

resVec   = zeros(iMax/iOut+1,1);
errorVec = zeros(iMax/iOut+1,1);

%%%%%%%
resVec(1)   = 1;
errorVec(1) = computeNormError2D_DG(mesh, dofm, x);
%%%%%%%

flag = 0;
i = 1;
while(i <= iMax)
    
    Ap = A*p;
    alpha = zzold/(Ap'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    z = A'*r;
    rrnew = r'*r;
    zznew = z'*z;
    
    %%%%%%%
    if(mod(i,iOut) == 0)
        resVec(i/iOut+1)   = sqrt(rrnew/rrini);
        errorVec(i/iOut+1) = computeNormError2D_DG(mesh, dofm, x);
        fprintf('[%i] %g %g\n', i, resVec(i/iOut+1), errorVec(i/iOut+1));
    end
    %%%%%%%
    
    if(sqrt(rrnew/rrini) < tol)
        flag = 1;
        return;
    end
    p = z + (zznew/zzold)*p;
    zzold = zznew;
    i = i+1;
end

end