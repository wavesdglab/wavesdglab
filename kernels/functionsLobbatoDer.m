% Assumption: x in [-1,1]

function val = functionsLobbatoDer(x,degree)

legendre = functionsLegendre(x,degree);
x = x(:);

N = degree+1;
val = zeros(size(x(:),1),N);

% nodal modes
val(:,1) = -1/2;
val(:,2) = +1/2;

% edge modes
for n=2:(N-1)
    val(:,n+1) = sqrt(n - 1/2) * legendre(:,n);
end

% N = degree+1;
% x = x(:);
% val = zeros(size(x,1),N);
% functionsKernelDer = functionsJacobiDer(x,1,1,N-2);
% 
% % nodal modes
% val(:,1) = -0.5;
% val(:,2) =  0.5;
% 
% % edge modes
% for n=1:N-2
%     val(:,n+2) = 0.25 * functionsKernelDer(:,n);
% end

end