function [souU, souDx, souDy, souVx, souVy] = mySourceSurface(x,y)

global M rho c omega;

souU = zeros(size(x));
souDx = zeros(size(x));
souDy = zeros(size(x));
souVx = zeros(size(x));
souVy = zeros(size(x));

% souU = exp(1i*k*x);

% souU = k^2*exp(1i*2*k*M/(M^2-1)*x);

% souU = exp(1i*k/(M+1)*x);

souU = exp(1i*omega/c*x);
% souVx = 1/(rho*c)*(1+(1-M)*(exp(1i*omega/c*x)-1));
souVx = 1/(rho*c)*exp(1i*omega/(c*M)*x)-M/(M-1)*exp(1i*omega/c*x)+M/(M-1)*exp(1i*omega/(c*M)*x);


A = (1+M)/2;
B = (1-M)/2;
souU = A*exp(1i*omega/(c*(1+M))*x)+B*exp(-1i*omega/(c*(1-M))*x);




% J = 8;
% N0 = floor(k*l/(pi*sqrt(1-M^2)));

% if abs(max(y)-min(y))<1e-6
%     souDy = zeros(size(x));
% 
% if abs(max(x)-min(x))<1e-6 && min(x)>1      % \Sigma_+
%     for j=1:J+1
%         if j==1
%             phi(j) = sqrt(1/l);
%         else
%             phi(j) = sqrt(2/l)*cos((n-1)*pi*y/l);
%         end
%         if (j<=N0)
%             beta(j) = (-k*M+sqrt(k^2-(n-1)^2*pi^2*(1-M)^2/l^2))/(1-M^2);
%         else
%             beta(j) = (-k*M+1i*sqrt((n-1)^2*pi^2*(1-M)^2/l^2-k^2))/(1-M^2);
%         end
%         souDx = souDx + 1i*beta(j)*phi(j);      % i still need to add (p,\phi_n)_{L^2(\Sigma_+)}
%     end
% end
% 
% if abs(max(x)-min(x))<1e-6 && min(x)<1      % \Sigma_-
%     for j=1:J+1
%         if j==1
%             phi(j) = sqrt(1/l);
%         else
%             phi(j) = sqrt(2/l)*cos((n-1)*pi*y/l);
%         end
%         if (j<=N0)
%             beta(j) = (-k*M-sqrt(k^2-(n-1)^2*pi^2*(1-M)^2/l^2))/(1-M^2);
%         else
%             beta(j) = (-k*M-1i*sqrt((n-1)^2*pi^2*(1-M)^2/l^2-k^2))/(1-M^2);
%         end
%         souDx = souDx + 1i*beta(j)*phi(j);      % i still need to add (p,\phi_n)_{L^2(\Sigma_-)}
%     end
% end



end