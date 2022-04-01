function [x,flag,relres,iter,resvec] = richardsonResPhy(A,b,tol,maxit,alpha,matA11,matA12,matA21,matA22,rhsA1,rhsA2)

% x = zeros(size(A,1),1);
% M = diag(diag(A));
% %M = sparse(1:size(A,1),1:size(A,1),1);

% normb = norm(b);
% iter = 0;
% while(iter < maxit)
%     iter = iter+1;
%     xNew = M\((M-A)*x + b);
%     x = alpha*xNew + (1-alpha)*x;
%     resvec = A*x - b;
%     relres = norm(A*x-b)/normb;
%     if(relres < tol)
%         flag = 1;
%         break;
%     end
% end
% flag = 2;

x = zeros(size(A,1),1);
%M = diag(diag(A));
M = sparse(1:size(A,1),1:size(A,1),1);

matA = matA11-(matA12/matA22)*matA21;
rhsA = rhsA1-(matA12/matA22)*rhsA2;

solA = matA11\(rhsA1 - matA12*x);
resA = matA*solA - rhsA;
rsiniA = resA'*resA;

normb = norm(b);
reshist = zeros(maxit,1);
iter = 0;
while(iter < maxit)
    iter = iter+1;
    xNew = M\((M-A)*x + b);
    x = alpha*xNew + (1-alpha)*x;
    resvec = A*x - b;
    relres = norm(A*x-b)/normb;
    
    solA = matA11\(rhsA1 - matA12*x);
    resA = matA*solA - rhsA;
    rsnewA = resA'*resA;
    
    reshist(iter) = sqrt(rsnewA/rsiniA);
    %disp(['                ' num2str(resvec(iter))]);
    if(reshist(iter) < tol)
        relres = reshist(iter);
        flag = 1;
        return;
    end
end
relres = reshist(maxit-1);
flag = 2;

end