function [resRedVec, resPhyVec, errorVec, i, flag] = solverCGNredu_CG(mesh, dofm, sys, tol, iMax, iOut)

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

%%%%%%%
solG = x;
solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
solA = [ solG ; solI ];
resRed = sys.rhsS - sys.matS*solG;
resPhy = sys.rhsA - sys.matA*solA;
resRedIni = resRed'*resRed;
resPhyIni = resPhy'*resPhy;
% resVec(1)    = 1;
resRedVec(1) = 1;
resPhyVec(1) = 1;
errorVec(1)  = computeNormError2D_CG(mesh, dofm, solA);
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
        solG = x;
        solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
        solA = [ solG ; solI ];
        resRed = sys.rhsS - sys.matS*solG;
        resPhy = sys.rhsA - sys.matA*solA;
        resRedNew = resRed'*resRed;
        resPhyNew = resPhy'*resPhy;
        % resVec(i/iOut+1)    = sqrt(rrnew/rrini);
        resRedVec(i/iOut+1) = sqrt(resRedNew/resRedIni);
        resPhyVec(i/iOut+1) = sqrt(resPhyNew/resPhyIni);
        errorVec(i/iOut+1)  = computeNormError2D_CG(mesh, dofm, solA);
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