% Assumption: x in [-1,1]

function val = functionsLagrangeLIN(x,degree)

x = x(:);
t = (1+x)/2;

N = degree+1;
val = zeros(size(x,1),N);

% P1
if(degree == 1)
    val(:,1) = 1-t;
    val(:,2) = t;
end

% P2
if(degree == 2)
    val(:,1) = (2*t-1).*(t-1);
    val(:,2) = t.*(2*t-1);
    val(:,3) = 4*t.*(1-t);
end

% P3
if(degree == 3)
    val(:,1) = -(1/2).*   (3*t-1).*(3*t-2).*(t-1);
    val(:,2) =  (9/2).*t         .*(3*t-2).*(t-1);
    val(:,3) = -(9/2).*t.*(3*t-1)         .*(t-1);
    val(:,4) =  (1/2).*t.*(3*t-1).*(3*t-2)       ;
end

if(degree > 3)
    error('Error: functionsLagrangeLIN() not available for degree P=%i', degree);
end

end