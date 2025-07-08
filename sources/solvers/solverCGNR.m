% Preconditined CGNR

function [resVec, errorVec, i, flag, x] = solverCGNR(mesh, dofm, sys, tol, iMax, iOut, computeError)

A = sys.matA;
b = sys.rhsA;
P = sys.matP;

% Right-preconditioning
x = zeros(size(A,2),1);
r = b - A*x;
z = P'\(A'*r);
p = z;
rrini = r'*r;
zzini = z'*z;
zzold = zzini;

% Symmetric-preconditioning
% x = zeros(size(A,2),1);
% r = b - A*x;
% z = P\(A'*(P\r));
% p = z;
% rrini = r'*r;
% zzini = z'*P*z;
% zzold = zzini;

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
    
    % Right-preconditioning
    q = A*(P\p);
    alpha = zzold/(q'*q);
    x = x + alpha*(P\p);
    r = r - alpha*q;
    z = P'\(A'*r);
    rrnew = r'*r;
    zznew = z'*z;
    beta = zznew/zzold;
    p = z + beta*p;
    zzold = zznew;
    
    % Symmetric-preconditioning
%     q = P\(A*p);
%     alpha = zzold/(q'*P*q);
%     x = x + alpha*p;
%     r = r - alpha*A*p;
%     z = P\(A'*(P\r));
%     rrnew = r'*r;
%     zznew = z'*P*z;
%     beta = zznew/zzold;
%     p = z + beta*p;
%     zzold = zznew;
    
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