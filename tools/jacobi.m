function [x,flag,relres,iter,resvec] = jacobi(A,b,tol,maxit,alpha)

x = zeros(size(A,1),1);
M = diag(diag(A));
%M = sparse(1:size(A,1),1:size(A,1),1);

normb = norm(b);
iter = 0;
while(iter < maxit)
    iter = iter+1;
    xNew = M\((M-A)*x + b);
    x = alpha*xNew + (1-alpha)*x;
    resvec = A*x - b;
    relres = norm(A*x-b)/normb;
    if(relres < tol)
        flag = 1;
        break;
    end
end
flag = 2;

end