function [x,flag,relres,iter,resvec] = solverCGN(A,b,resTol,maxit)


%try chol(A'*A);
%    disp('Matrix is symmetric positive definite.')
%catch ME
%    disp('Matrix is not symmetric positive definite')
%end

% [x,flag,relres,iter,resvec] = solverCG(A'*A,A'*b,resTol,maxit);
% size(A)
% maxit
[x,flag,relres,iter,resvec] = pcg(A'*A,A'*b,resTol,maxit);


% x = zeros(size(A,1),1);
% r = A'*b - A'*A*x;
% p = r;
% rsini = r'*r;
% rsold = rsini;
% 
% resvec = zeros(maxit,1);
% iter = 1;
% while(iter <= maxit)
%     Ap = A'*A*p;
%     alpha = rsold/(p'*Ap);
%     x = x + alpha*p;
%     r = r - alpha*Ap;
%     rsnew = r'*r;
%     resvec(iter) = sqrt(rsnew/rsini);
%     %disp(['                ' num2str(resvec(iter))]);
%     if(resvec(iter) < resTol)
%         relres = resvec(iter);
%         flag = 1;
%         return;
%     end
%     p = r + (rsnew/rsold)*p;
%     rsold = rsnew;
%     iter = iter+1;
% end
% relres = resvec(maxit-1);
% flag = 2;



%disp(['    Min(eigval A^TA)   ' num2str(min(eigs(A'*A)))]);
%disp(['    Max(eigval A^TA)   ' num2str(max(eigs(A'*A)))]);
%disp(['    Cond(A^TA)         ' num2str(condest(A'*A))]);

end
