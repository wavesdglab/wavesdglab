% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

% Richardson with symmetric preconditioning

function [resRedVec] = Richardson(sys, tol, iMax, iOut, alpha)

A = sys.matS;
b = sys.rhsS;
Pinv = sys.matPinv;

x = zeros(size(A,2),1);
r = b-A*x;
rrini = r'*r;

resRedVec = zeros(iMax/iOut+1,1); 

%%%%%
resRedVec(1) = 1;
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
        resRedVec(i/iOut+1) = sqrt(rrnew/rrini);
        fprintf('[%i] %g\n', i, resRedVec(i/iOut+1));
    end
    %%%%%%%
    
    if(sqrt(rrnew/rrini) < tol)
        flag = 1;
        return;
    end

    i = i+1;
end

end