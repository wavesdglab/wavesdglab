% CGNR iteration with symmetric preconditioning

function [resRedVec, resPhyVec, errorVec, i, flag, xPhy] = solverCGNredu_DG(mesh, dofm, sys, tol, iMax, iOut, computeError)

A = sys.matS;
b = sys.rhsS;
P = sys.matP;
Pinv = sys.matPinv;

x = zeros(size(A,2),1);
r = b-A*x;
s = Pinv*r;
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
xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
rPhy = sys.rhsPhy - sys.matPhy*xPhy;
resPhyIni = rPhy'*rPhy;
resRedVec(1) = 1;
resPhyVec(1) = 1;
errorVec(1) = computeError(mesh, dofm, xPhy);
%%%%%%%

flag = 0;
i = 1;
while(i <= iMax)
    
%     if(mod(i,iOut) == 0)
%         [x,flag] = pcg(A'*A,A'*b,tol,i);
%         r = b-A*x;
%         xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
%         rPhy = sys.rhsPhy - sys.matPhy*xPhy;
%         resPhyNew = rPhy'*rPhy;
%         resRedVec(i/iOut+1) = sqrt(r'*r/rrini);
%         resPhyVec(i/iOut+1) = sqrt(resPhyNew/resPhyIni);
%         errorVec(i/iOut+1) = computeNormError2D_DG(mesh, dofm, xPhy);
%         fprintf('[%i] %g %g\n', i, resRedVec(i/iOut+1), errorVec(i/iOut+1));
%     end
    
    v = A*p;
    w = Pinv*v;
    alpha = zzold/(v'*w);
    x = x + alpha*p;
    r = r - alpha*v;
    s = Pinv*r;
    y = A'*s;
    z = Pinv*y;
    rrnew = r'*r;
    zznew = y'*z;
    p = z + (zznew/zzold)*p;
    zzold = zznew;
    
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

end