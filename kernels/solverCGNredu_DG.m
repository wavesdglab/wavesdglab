function [resRedVec, resPhyVec, errorVec, i, flag] = solverCGNredu_DG(mesh, dofm, sys, tol, iMax, iOut)

A = sys.matS;
b = sys.rhsS;
x = zeros(size(A,2),1);
r = b - A*x;

z = A'*r;
p = z;
rrini = r'*r;
zzini = z'*z;
zzold = zzini;

% resVec       = zeros(iMax/iOut+1,1);
resRedVec    = zeros(iMax/iOut+1,1);
resPhyVec    = zeros(iMax/iOut+1,1);
errorVec     = zeros(iMax/iOut+1,1);
%errorPostVec = zeros(iMax/iOut+1,1);

%%%%%%%
solG = sys.precR*x;
solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
%[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solI);
resPhy = sys.rhsPhy - sys.matPhy*solI;
resPhyIni = resPhy'*resPhy;
% resVec(1)       = 1;
resRedVec(1)    = 1;
resPhyVec(1)    = 1;
errorVec(1)     = computeNormError2D_DG(mesh, dofm, solI);
%errorPostVec(1) = computeNormError2D_DG(mesh, dofmPost, solApost);
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
        solG = sys.precR*x;
        solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
        %[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solI);
        resPhy = sys.rhsPhy - sys.matPhy*solI;
        resPhyNew = resPhy'*resPhy;
        % resVec(i/iOut+1)       = sqrt(zznew/zzini);
        resRedVec(i/iOut+1)    = sqrt(rrnew/rrini);
        resPhyVec(i/iOut+1)    = sqrt(resPhyNew/resPhyIni);
        errorVec(i/iOut+1)     = computeNormError2D_DG(mesh, dofm, solI);
        %errorPostVec(i/iOut+1) = computeNormError2D_DG(mesh, dofmPost, solApost);
        fprintf('[%i] %g %g\n', i, resRedVec(i/iOut+1), errorVec(i/iOut+1));
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