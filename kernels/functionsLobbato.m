% Assumption: x in [-1,1]

function val = functionsLobbato(x,N)

% val = functionsBernstein(x,N);

legendreInt = functionsLegendreInt(x,N);
x = x(:);

val = zeros(size(x(:),1),N);
val(:,1) = (1-x)/2;  % order 0
val(:,2) = (1+x)/2;  % order 1
for n=2:(N-1)
    val(:,n+1) = sqrt(n - 1/2) * legendreInt(:,n);  % order n
end

end