function [resVec, errorVec, iter, flag] = solverCGN_DG(mesh, dofm, sys, tol, iMax, iOut)

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
iter = iMax;
for i=1:iMax
    
%     if(mod(i,iOut) == 0)
%         [x,flag] = pcg(A'*A,A'*b,tol,i);
%         r = b - A*x;
%         resVec(i/iOut+1) = sqrt(r'*r/rrini);
%         errorVec(i/iOut+1) = computeNormError2D_DG(mesh, dofm, x);
%         fprintf('[%i] %g %g\n', i, resVec(i/iOut+1), errorVec(i/iOut+1));
%     end
    
    Ap = A*p;
    alpha = zzold/(Ap'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    z = A'*r;
    rrnew = r'*r;
    zznew = z'*z;
    relRes = sqrt(rrnew/rrini);
    
    %%%%%%%
    if(mod(i,iOut) == 0)
        resVec(i/iOut+1) = relRes;
        errorVec(i/iOut+1) = computeNormError2D_DG(mesh, dofm, x);
        fprintf('[%i] %g %g\n', i, resVec(i/iOut+1), errorVec(i/iOut+1));
    end
    %%%%%%%
    
    p = z + (zznew/zzold)*p;
    zzold = zznew;
    
    if (relRes <= tol)
        iter = i;
        flag = 1;
        break;
    end
end
end