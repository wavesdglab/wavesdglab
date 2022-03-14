function [x,flag,relres,iter,resvec] = conjgradn(A,b,resTol,maxit)


%try chol(A'*A);
%    disp('Matrix is symmetric positive definite.')
%catch ME
%    disp('Matrix is not symmetric positive definite')
%end

%[x,flag,relres,iter,resvec] = conjgrad(A'*A,A'*b,resTol,maxit);
[x,flag,relres,iter,resvec] = pcg(A'*A,A'*b,resTol,maxit);

%disp(['    Min(eigval A^TA)   ' num2str(min(eigs(A'*A)))]);
%disp(['    Max(eigval A^TA)   ' num2str(max(eigs(A'*A)))]);
%disp(['    Cond(A^TA)         ' num2str(condest(A'*A))]);

end
