function [resVec, resRedVec, resPhyVec, error, iter, flag] = solverCGNredu_CG(mesh, dofm, sys, tol, maxit, itoutput)

A = sys.matS'*sys.matS;
b = sys.matS'*sys.rhsS;

x = zeros(size(A,1),1);
r = b - A*x;
p = r;
rsini = r'*r;
rsold = rsini;

%%%%%
solG = x;
solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
solA = [ solG ; solI ];
resRed = sys.rhsS - sys.matS*solG;
resPhy = sys.rhsA - sys.matA*solA;
resRedIni = resRed'*resRed;
resPhyIni = resPhy'*resPhy;
%%%%%

resVec  = zeros(maxit/itoutput+1,1);
resRedVec = zeros(maxit/itoutput+1,1);
resPhyVec = zeros(maxit/itoutput+1,1);
error   = zeros(maxit/itoutput+1,1);

resVec(1) = 1;
resRedVec(1) = 1;
resPhyVec(1) = 1;
error(1) = computeNormError2D_CG(mesh, dofm, solA);

iter = 1;
while(iter <= maxit)
    Ap = A*p;
    alpha = rsold/(p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    rsnew = r'*r;
    resVecNew = sqrt(rsnew/rsini);
    
    %%%%%%%
    if(mod(iter,itoutput)==0)
        solG = x;
        solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
        solA = [ solG ; solI ];
        resRed = sys.rhsS - sys.matS*solG;
        resPhy = sys.rhsA - sys.matA*solA;
        resRedNew = resRed'*resRed;
        resPhyNew = resPhy'*resPhy;
        resVec(iter/itoutput+1) = resVecNew;
        resRedVec(iter/itoutput+1) = sqrt(resRedNew/resRedIni);
        resPhyVec(iter/itoutput+1) = sqrt(resPhyNew/resPhyIni);
        error(iter/itoutput+1) = computeNormError2D_CG(mesh, dofm, solA);
    end
    %%%%%%%
    
    %disp(['                ' num2str(resvec(iter))]);
    if(resVecNew < tol)
        flag = 1;
        return;
    end
    p = r + (rsnew/rsold)*p;
    rsold = rsnew;
    iter = iter+1;
end
flag = 2;

end
