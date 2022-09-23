function [resRedVec, resPhyVec, errorVec, i, flag] = solverRichardson_DG(mesh, dofm, sys, tol, iMax, iOut, alpha)

A = sys.matS;
b = sys.rhsS;
x = zeros(size(A,2),1);
r = b - A*x;

rrini = r'*r;

resRedVec    = zeros(iMax/iOut+1,1); 
resPhyVec    = zeros(iMax/iOut+1,1); 
errorVec     = zeros(iMax/iOut+1,1);
%errorPostVec = zeros(iMax/iOut+1,1);

%%%%%
solG = sys.precR*x;
solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
%[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solI);
resPhy = sys.rhsPhy - sys.matPhy*solI;
resPhyIni = resPhy'*resPhy;
resRedVec(1)    = 1;
resPhyVec(1)    = 1;
errorVec(1)     = computeNormError2D_DG(mesh, dofm, solI);
%errorPostVec(1) = computeNormError2D_DG(mesh, dofmPost, solApost);
%%%%%

flag = 0;
i = 1;
while(i <= iMax)
    
    xNew = x - A*x + b;
    x = alpha*xNew + (1-alpha)*x;
    r = b - A*x;
    rrnew = r'*r;
    
    %%%%%%%
    if(mod(i,iOut)==0)
        solG = sys.precR*x;
        solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
        %[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solI);
        resPhy = sys.rhsPhy - sys.matPhy*solI;
        resPhyNew = resPhy'*resPhy;
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
    
    i = i+1;
end

end