function [x,flag,relres,iter,resvec] = conjgrad(A,b,tol,maxit)

%minEig = min(eigs(A))
%cond = condest(A)

x = zeros(size(A,1),1);
r = b - A*x;
p = r;
rsini = r'*r;
rsold = rsini;

resvec = zeros(maxit,1);
iter = 1;
while(iter <= maxit)
    Ap = A*p;
    alpha = rsold/(p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    rsnew = r'*r;
    resvec(iter) = sqrt(rsnew/rsini);
    %disp(['                ' num2str(resvec(iter))]);
    if(resvec(iter) < tol)
        relres = resvec(iter);
        flag = 1;
        return;
    end
    p = r + (rsnew/rsold)*p;
    rsold = rsnew;
    iter = iter+1;
end
relres = resvec(maxit-1);
flag = 2;

end
