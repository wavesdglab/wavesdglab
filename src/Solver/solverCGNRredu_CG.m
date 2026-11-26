% CGNR with right preconditioning

function [resRedVec, resPhyVec, errorVec, i, flag, xPhy] = solverCGNRredu_CG(mesh, dofm, sys, tol, iMax, iOut, computeError)

A = sys.matS;
b = sys.rhsS;
P = 1;

% Right-preconditioning
x = zeros(size(A,2),1);
r = b - A*x;
z = P'\(A'*r);
p = z;
rrini = r'*r;
zzini = z'*z;
zzold = zzini;

resRedVec = zeros(iMax/iOut+1,1);
resPhyVec = zeros(iMax/iOut+1,1);
errorVec  = zeros(iMax/iOut+1,1);

%%%%%%%
xPhy = [ x ; sys.matII\(sys.rhsI-sys.matIG*x) ];
rPhy = sys.rhsA - sys.matA*xPhy;
resPhyIni = sqrt(rPhy'*rPhy);
resRedVec(1) = 1;
resPhyVec(1) = 1;
errorVec(1) = computeError(mesh, dofm, xPhy);
fprintf('[%i] %g %g %g\n', 1, resRedVec(1), resPhyVec(1), errorVec(1));
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
    
    %%%%%%%
    if(mod(i,iOut) == 0)
        xPhy = [ x ; sys.matII\(sys.rhsI-sys.matIG*x) ];
        rPhy = sys.rhsA - sys.matA*xPhy;
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

xPhy = [ x ; sys.matII\(sys.rhsI-sys.matIG*x) ];

end