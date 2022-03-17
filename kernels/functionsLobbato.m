% Assumption: x in [-1,1]

function val = functionsLobbato(x,degree)
                                                    
% legendreInt = functionsLegendreInt(x,degree);
% x = x(:);
% 
% N = degree+1;
% val = zeros(size(x(:),1),N);
% 
% % nodal modes
% val(:,1) = (1-x)/2;
% val(:,2) = (1+x)/2;
% 
% % edge modes
% for n=2:(N-1)
%     val(:,n+1) = sqrt(n - 1/2) * legendreInt(:,n);
% end

N = degree+1;
x = x(:);
val = zeros(size(x,1),N);
functionsKernel = functionsJacobi(x,1,1,N-2);

% nodal modes
val(:,1) = (1-x)/2;
val(:,2) = (1+x)/2;

% edge modes
for n=1:N-2
    val(:,n+2) = val(:,1) .* val(:,2) .* functionsKernel(:,n);
end

end