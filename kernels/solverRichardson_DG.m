function [resRedVec, resPhyVec, error, errorPP, iter, flag] = solverRichardson_DG(mesh, dofm, sys, tol, maxit, itoutput, alpha)

A = sys.matS;
b = sys.rhsS;
x = zeros(size(A,1),1);

%%%%%
solG = x;
solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
resRed = sys.rhsS - sys.matS*solG;
resPhy = sys.rhsPhy - sys.matPhy*solI;
resRedIni = resRed'*resRed;
resPhyIni = resPhy'*resPhy;
%%%%%

resRedVec = zeros(maxit/itoutput+1,1); 
resPhyVec = zeros(maxit/itoutput+1,1); 
error   = zeros(maxit/itoutput+1,1);
errorPP = zeros(maxit/itoutput+1,1);

resRedVec(1) = 1;
resPhyVec(1) = 1;
error(1) = computeNormError2D_DG(mesh, dofm, solI);
[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solI);
errorPP(1) = computeNormError2D_DG(mesh, dofmPost, solApost);

iter = 1;
while(iter < maxit+1)
    if(mod(iter/itoutput,1)==0)
        disp(num2str(iter))
    end
    
    xNew = x - A*x + b;
    x = alpha*xNew + (1-alpha)*x;
    resRed = b - A*x;
    resRedNew = resRed'*resRed;
    
    %%%%%%%
    if(mod(iter,itoutput)==0)
        solG = x;
        solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
        resPhy = sys.rhsPhy - sys.matPhy*solI;
        resPhyNew = resPhy'*resPhy;
        resRedVec(iter/itoutput+1) = sqrt(resRedNew/resRedIni);
        resPhyVec(iter/itoutput+1) = sqrt(resPhyNew/resPhyIni);
        error(iter/itoutput+1) = computeNormError2D_DG(mesh, dofm, solI);
        [solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solI);
        errorPP(iter/itoutput+1) = computeNormError2D_DG(mesh, dofmPost, solApost);
    end
    %%%%%%%
    
    if(sqrt(resRedNew/resRedIni) < tol)
        flag = 1;
        return;
    end
    
    iter = iter+1;
end
flag = 2;

end