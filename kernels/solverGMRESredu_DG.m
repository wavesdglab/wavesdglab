function [resRedVec, resPhyVec, errorVec, i, flag] = solverGMRESredu_DG(mesh, dofm, sys, tol, iMax, iOut)

A = sys.matS;
b = sys.rhsS;
x = zeros(size(A,2),1);
r = b - A*x;

sn = zeros(iMax,1);
cs = zeros(iMax,1);
e1 = zeros(iMax+1,1);
e1(1) = 1;
r_norm = norm(r);
H = zeros(iMax+1,iMax+1);
Q = zeros(size(A,2),iMax+1);
Q(:,1) = r / r_norm;
beta = r_norm * e1;

resRedVec    = zeros(iMax/iOut+1,1);
resPhyVec    = zeros(iMax/iOut+1,1);
errorVec     = zeros(iMax/iOut+1,1);

%%%%%%%
solG = sys.precR*x;
solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
%[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solI);
resPhy = sys.rhsPhy - sys.matPhy*solI;
resPhyIni = resPhy'*resPhy;
resRedVec(1)    = 1;
resPhyVec(1)    = 1;
errorVec(1)     = computeNormError2D_DG(mesh, dofm, solI);
%errorPostVec(1) = computeNormError2D_DG(mesh, dofmPost, solApost);
%%%%%%%

flag = 0;
i = 1;
while(i <= iMax)
    
    % Arnoldi iteration – Add one vector to basis Q and orthogonalize it
    Q(:,i+1) = A*Q(:,i);
    for j = 1:i
        H(j,i) = Q(:,i+1)' * Q(:,j);
        Q(:,i+1) = Q(:,i+1) - H(j,i) * Q(:,j);
    end
    H(i+1,i) = norm(Q(:,i+1));
    Q(:,i+1) = Q(:,i+1) / H(i+1,i);
    
    % Apply for ith column
    for j = 1:i-1
        tmp      =  cs(j) * H(j,i) + sn(j) * H(j+1,i);
        H(j+1,i) = -sn(j) * H(j,i) + cs(j) * H(j+1,i);
        H(j,i)   = tmp;
    end
    
    % Update the next sin cos values for rotation
    tmp = sqrt(H(i,i)^2 + H(i+1,i)^2);
    cs(i) = H(i,i)/tmp;
    sn(i) = H(i+1,i)/tmp;
    
    % Eliminate H(j+1,j)
    H(i,i) = cs(i) * H(i,i) + sn(i) * H(i+1,i);
    H(i+1,i) = 0;
    
    % Update the residual vector
    beta(i+1) = -sn(i) * beta(i);
    beta(i)   = cs(i) * beta(i);
    error     = abs(beta(i+1)) / norm(b);
    
    %%%%%%%
    if(mod(i,iOut) == 0)
        
        y = H(1:i,1:i) \ beta(1:i);
        x = Q(:,1:i) * y;
        
        solG = sys.precR*x;
        solI = sys.matIIinv*(sys.rhsI-sys.matIG*solG);
        %[solApost, dofmPost] = computeSolPostPro2D_DG(mesh, dofm, solI);
        resPhy = sys.rhsPhy - sys.matPhy*solI;
        resPhyNew = resPhy'*resPhy;
        resRedVec(i/iOut+1)    = error;
        resPhyVec(i/iOut+1)    = sqrt(resPhyNew/resPhyIni);
        errorVec(i/iOut+1)     = computeNormError2D_DG(mesh, dofm, solI);
        %errorPostVec(i/iOut+1) = computeNormError2D_DG(mesh, dofmPost, solApost);
        fprintf('[%i] %g %g\n', i, resRedVec(i/iOut+1), errorVec(i/iOut+1));
    end
    %%%%%%%
    
    if (error <= tol)
        flag = 1;
        break;
    end
    
    i = i+1;
end

end