function [resVec, resRedVec, resPhyVec, error, errorPP, iter, flag] = solverCGNreduDG(mesh, dofm, sys, tol, maxit, itoutput)

A = sys.matS;
b = sys.rhsS;

x = zeros(size(A,1),1);
r = b - A*x;
z = A'*r;
p = z;
rrini = r'*r;
zzini = z'*z;
zzold = zzini;

%%%%%
solG = x;
solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
resRed = sys.rhsS - sys.matS*solG;
resPhy = sys.rhsPhy - sys.matPhy*solI;
resRedIni = resRed'*resRed;
resPhyIni = resPhy'*resPhy;
%%%%%

resVec  = zeros(maxit/itoutput+1,1);
resRedVec = zeros(maxit/itoutput+1,1);
resPhyVec = zeros(maxit/itoutput+1,1);
error   = zeros(maxit/itoutput+1,1);
errorPP = zeros(maxit/itoutput+1,1);

resVec(1) = 1;
resRedVec(1) = 1;
resPhyVec(1) = 1;
error(1) = computeNormError2D_DG(mesh, dofm, solI);
[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solI);
errorPP(1) = computeNormError2D_DG(mesh, dofmPost, solApost);

iter = 1;

while(iter <= maxit)
    if(mod(iter,100)==0)
        disp(num2str(iter))
    end
    
    Ap = A*p;
    alpha = zzold/(Ap'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    z = A'*r;
    rrnew = r'*r;
    zznew = z'*z;
    
    %%%%%%%
    if(mod(iter,itoutput)==0)
        solG = x;
        solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
        resPhy = sys.rhsPhy - sys.matPhy*solI;
        resPhyNew = resPhy'*resPhy;
        resVec(iter/itoutput+1)    = sqrt(zznew/zzini);
        resRedVec(iter/itoutput+1) = sqrt(rrnew/rrini);
        resPhyVec(iter/itoutput+1) = sqrt(resPhyNew/resPhyIni);
        error(iter/itoutput+1) = computeNormError2D_DG(mesh, dofm, solI);
        [solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solI);
        errorPP(iter/itoutput+1) = computeNormError2D_DG(mesh, dofmPost, solApost);
    end
    %%%%%%%
    
    if(sqrt(rrnew/rrini) < tol)
        flag = 1;
        return;
    end
    
    p = z + (zznew/zzold)*p;
    zzold = zznew;
    iter = iter+1;
end
flag = 2;

end