function [x,flag,relRes,iter] = solverGMRES(A,b,tol,iMax)

[x,flag,relRes,iter] = gmres(A,b,[],tol,iMax);
MATL = [flag,relRes,iter(2)]

% ==========

x = zeros(size(A,2),1);
r = b - A*x;

sn = zeros(iMax,1);
cs = zeros(iMax,1);
H = zeros(iMax+1,iMax);
Q = zeros(size(A,2),iMax+1);
Q(:,1) = r/norm(r);
beta = zeros(iMax+1,1);
beta(1) = norm(r);

flag = 0;
iter = iMax;
relRes = 1;
for i = 1:iMax
    
    % Arnoldi iteration – Add one vector to basis Q and orthogonalize it
    Q(:,i+1) = A*Q(:,i);
    for j = 1:i
        H(j,i) = Q(:,j)' * Q(:,i+1);
        Q(:,i+1) = Q(:,i+1) - H(j,i)' * Q(:,j);
    end
    H(i+1,i) = norm(Q(:,i+1));
    Q(:,i+1) = Q(:,i+1) / H(i+1,i);
    
    % Apply previous Givens matrix to ith column
    for j = 1:i-1
        matGivens = [ cs(j)' sn(j)' ; -sn(j) cs(j) ];
        H(j:j+1,i) = matGivens * H(j:j+1,i);
    end
    
    % Compute the new Givens matrix
    tmp = sqrt(abs(H(i,i))^2 + H(i+1,i)^2);
    cs(i) = H(i,i)/tmp;    % complex
    sn(i) = H(i+1,i)/tmp;  % real
    matGivens = [ cs(i)' sn(i)' ; -sn(i) cs(i) ];
    
    % Apply new Givens matrix to ith column of H and residual vector
    H(i:i+1,i)  = matGivens * H(i:i+1,i);
    beta(i:i+1) = matGivens * beta(i:i+1);
    
    % Update the residual vector
    relRes = abs(beta(i+1)) / norm(beta(1));
    
    % Update the solution
    y = H(1:i,1:i) \ beta(1:i);
    x = Q(:,1:i) * y;
    
    if (relRes <= tol)
        iter = i;
        flag = 1;
        break;
    end
end

MINE = [flag,relRes,iter]

end