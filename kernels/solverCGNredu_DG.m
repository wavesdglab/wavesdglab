function [resRedVec, resPhyVec, errorVec, iter, flag] = solverCGNredu_DG(mesh, dofm, sys, tol, iMax, iOut)

A = sys.matS;
b = sys.rhsS;
x = zeros(size(A,2),1);
r = b - A*x;

z = A'*r;
p = z;
rrini = r'*r;
zzini = z'*z;
zzold = zzini;

resRedVec = zeros(iMax/iOut+1,1);
resPhyVec = zeros(iMax/iOut+1,1);
errorVec  = zeros(iMax/iOut+1,1);

%%%%%%%
solG = sys.precR*x;
solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
resPhy = sys.rhsPhy - sys.matPhy*solI;
resPhyIni = resPhy'*resPhy;
resRedVec(1) = 1;
resPhyVec(1) = 1;
errorVec(1) = computeNormError2D_DG(mesh, dofm, solI);
%%%%%%%

flag = 0;
iter = iMax;
for i=1:iMax
    
    if(mod(i,iOut) == 0)
        %[x,flag] = pcg(A'*A,A'*b,tol,i);
        [x,flag] = qmr(A,b,tol,i);
        %[x,flag] = qmr(A'*A,A'*b,tol,i);
        r = b - A*x;
        solG = sys.precR*x;
        solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
        resPhy = sys.rhsPhy - sys.matPhy*solI;
        resPhyNew = resPhy'*resPhy;
        resRedVec(i/iOut+1) = sqrt(r'*r/rrini);
        resPhyVec(i/iOut+1) = sqrt(resPhyNew/resPhyIni);
        errorVec(i/iOut+1) = computeNormError2D_DG(mesh, dofm, solI);
        fprintf('[%i] %g %g\n', i, resRedVec(i/iOut+1), errorVec(i/iOut+1));
    end
    
%     Ap = A*p;
%     alpha = zzold/(Ap'*Ap);
%     x = x + alpha*p;
%     r = r - alpha*Ap;
%     z = A'*r;
%     rrnew = r'*r;
%     zznew = z'*z;
%     relRes = sqrt(rrnew/rrini);
%     
%     %%%%%%%
%     if(mod(i,iOut) == 0)
%         solG = sys.precR*x;
%         solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
%         resPhy = sys.rhsPhy - sys.matPhy*solI;
%         resPhyNew = resPhy'*resPhy;
%         resRedVec(i/iOut+1) = relRes;
%         resPhyVec(i/iOut+1) = sqrt(resPhyNew/resPhyIni);
%         errorVec(i/iOut+1) = computeNormError2D_DG(mesh, dofm, solI);
%         fprintf('[%i] %g %g\n', i, resRedVec(i/iOut+1), errorVec(i/iOut+1));
%     end
%     %%%%%%%
%     
%     p = z + (zznew/zzold)*p;
%     zzold = zznew;
%     
%     if (relRes <= tol)
%         iter = i;
%         flag = 1;
%         break;
%     end
end
end