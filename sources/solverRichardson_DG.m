% Richardson iteration with symmetric preconditioning

function [resRedVec, resPhyVec, errorVec, i, flag, xPhy] = solverRichardson_DG(mesh, dofm, sys, tol, iMax, iOut, alpha, computeError)

A = sys.matS;
b = sys.rhsS;
P = sys.matP;
Pinv = sys.matPinv;
x = zeros(size(A,2),1);
r = b-A*x;

rrini = r'*r;

resRedVec = zeros(iMax/iOut+1,1); 
resPhyVec = zeros(iMax/iOut+1,1); 
errorVec  = zeros(iMax/iOut+1,1);

%%%%%
xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
rPhy = sys.rhsPhy - sys.matPhy*xPhy;
resPhyIni = rPhy'*rPhy;
resRedVec(1) = 1;
resPhyVec(1) = 1;
errorVec(1) = computeError(mesh, dofm, xPhy);
%%%%%

flag = 0;
i = 1;
while(i <= iMax)
    
    xNew = Pinv * (P*x - A*x + b);
    x = alpha*xNew + (1-alpha)*x;
    r = b - A*x;
    rrnew = r'*r;
    
    %%%%%%%
    if(mod(i,iOut)==0)
        xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
        rPhy = sys.rhsPhy - sys.matPhy*xPhy;
        resPhyNew = rPhy'*rPhy;
        resRedVec(i/iOut+1) = sqrt(rrnew/rrini);
        resPhyVec(i/iOut+1) = sqrt(resPhyNew/resPhyIni);
        errorVec(i/iOut+1) = computeError(mesh, dofm, xPhy);
        fprintf('[%i] %g %g\n', i, resRedVec(i/iOut+1), errorVec(i/iOut+1));
    end
    %%%%%%%
    
    if(sqrt(rrnew/rrini) < tol)
        flag = 1;
        xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);
        return;
    end
    
    i = i+1;
end

xPhy = sys.matIIinv*(sys.rhsI-sys.matIG*x);

end