% Assumption: x in [-1,1]

function val = functionsLagrangeDerLIN(x,degree)

x = x(:);
t = (1+x)/2;
tDX = 1/2;

N = degree+1;
val = zeros(size(x,1),N);

% P1
if(degree == 1)
    val(:,1) = -tDer;
    val(:,2) = tDer;
end

% P2
if(degree == 2)
    val(:,1) = 2*tDX.*(t-1) + (2*t-1).*tDX;
    val(:,2) = tDX.*(2*t-1) + t.*(2*tDX);
    val(:,3) = 4*tDX.*(1-t) + 4*t.*(-tDX);
end

% P3
if(degree == 3)
    val(:,1) = - (1/2).*(3*tDX).*(3*t-2).*(t-1) - (1/2).*(3*t-1).*(3*tDX).*(t-1) - (1/2).*(3*t-1).*(3*t-2).*(tDX);
    val(:,2) = + (9/2).*tDX.*(3*t-2).*(t-1) + (9/2).*t.*(3*tDX).*(t-1) + (9/2).*t.*(3*t-2).*(tDX);
    val(:,3) = - (9/2).*tDX.*(3*t-1).*(t-1) - (9/2).*t.*(3*tDX).*(t-1) - (9/2).*t.*(3*t-1).*(tDX);
    val(:,4) = + (1/2).*tDX.*(3*t-1).*(3*t-2) + (1/2).*t.*(3*tDX).*(3*t-2) + (1/2).*t.*(3*t-1).*(3*tDX);
end

end