function [x,flag,relres,iter,resvec] = conjgradnResPhy(A,b,resTol,maxit,matA11,matA12,matA21,matA22,rhsA1,rhsA2)

x = zeros(size(A,1),1);
r = A'*b - A'*A*x;
p = r;
rsini = r'*r;
rsold = rsini;

matA = matA11-(matA12/matA22)*matA21;
rhsA = rhsA1-(matA12/matA22)*rhsA2;

solA = matA11\(rhsA1 - matA12*x);
resA = matA*solA - rhsA;
rsiniA = resA'*resA;

resvec = zeros(maxit,1);
iter = 1;
while(iter <= maxit)
    Ap = A'*A*p;
    alpha = rsold/(p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    rsnew = r'*r;
    
    solA = matA11\(rhsA1 - matA12*x);
    resA = matA*solA - rhsA;
    rsnewA = resA'*resA;
    
    resvec(iter) = sqrt(rsnewA/rsiniA);
    %disp(['                ' num2str(resvec(iter))]);
    if(resvec(iter) < resTol)
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
