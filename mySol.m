function [sol, solDx, solDy, solF] = mySol(x,y)

global k;

% sol   = exp(1i*k*x);
% solF  = 0*x;
% solDx = 1i*omega * sol;
% solDy = 0*x;

% sol   = -1/k * cos(k*(1-x)) / sin(2*k);
% solF  = 0*x;  % [- Delta - k^2] u
% solDx = -sin(k*(1-x)) / sin(2*k);
% solDy = 0*x;

% sol   = sin(pi*x/2) .* sin(pi*y/2);
% solF  = (0.5*pi^2 - k^2) * sol;  % [- Delta - omega^2] u
% solDx = 0*x;
% solDy = 0*x;

% Dirichlet problem
% m = 3;
% n = 2;
% sol   = sin(m*pi*x) .* sin(n*pi*y);
% solF  = ((m*m+n*n)*pi^2 - k^2) * sol;
% solDx = m*pi * cos(m*pi*x) .* sin(n*pi*y);
% solDy = n*pi * sin(m*pi*x) .* cos(n*pi*y);

% Neumann problem
m = 2;
n = 3;
sol   = cos(m*pi*x) .* cos(n*pi*y);
solF  = ((m*m+n*n)*pi^2 - k^2) * sol;
solDx = -m*pi * sin(m*pi*x) .* cos(n*pi*y);
solDy = -n*pi * cos(m*pi*x) .* sin(n*pi*y);

% Duct problem
% m = 6;
% kx = sqrt(k^2 - (m*pi)^2);
% alpha = (k-kx)/(k+kx) * exp(-4i*kx);
% sol   = 1/(1i*kx) * cos(m*pi*y) .* (exp(-1i*kx*x) - alpha * exp(1i*kx*x))/(1+alpha);
% solF  = 0*x;
% solDx = 1/(1i*kx) * cos(m*pi*y) .* (-1i*kx*exp(-1i*kx*x) - 1i*kx*alpha * exp(1i*kx*x))/(1+alpha);
% solDy = - m*pi * 1/(1i*kx) * sin(m*pi*y) .* (exp(-1i*kx*x) - alpha * exp(1i*kx*x))/(1+alpha);

end