% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Richardson with symmetric preconditioning

function [resRedVec, resPhyVec, errorVec, i, flag, x] = solverRichardson_DG(mesh, dofm, sys, tol, iMax, iOut, alpha, computeError)

A = sys.matA;
b = sys.rhsA;
Pinv = sys.matPinv;

x = zeros(size(A,2),1);
r = b-A*x;
rrini = r'*r;

resRedVec = zeros(iMax/iOut+1,1); 
resPhyVec = zeros(iMax/iOut+1,1); 
errorVec  = zeros(iMax/iOut+1,1);

%%%%%
resPhyIni = r'*r;
resRedVec(1) = 1;
resPhyVec(1) = 1;
errorVec(1) = computeError(mesh, dofm, x);
%%%%%

flag = 0;
i = 1;
while(i <= iMax)
    
    % xNew = Pinv * (P*x - A*x + b);
    % x = alpha*xNew + (1-alpha)*x;
    
    x = alpha*Pinv*r + x;
    r = b-A*x;
    rrnew = r'*r;
    
    %%%%%%%
    if(mod(i,iOut)==0)
        resPhyNew = r'*r;
        resRedVec(i/iOut+1) = sqrt(rrnew/rrini);
        resPhyVec(i/iOut+1) = sqrt(resPhyNew/resPhyIni);
        errorVec(i/iOut+1) = computeError(mesh, dofm, x);
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