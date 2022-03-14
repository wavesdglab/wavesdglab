function output = fovals(A,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Wilmer Henao    wi-henao@uniandes.edu.co
%   Department of Mathematics
%   Universidad de los Andes
%   Colombia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

EPS = 0.0000000000000000000001*i;
[m,n] = size(A);
if m ~= n
   disp('The matrix must be square')
   return
end

tr = trace(A)/m;
A = A - tr.*eye(m);

if nargin == 1
    k = 500;
end

pto = zeros(1,k);
for ind = 1:k
    theta = 2*pi*ind/k;
    [vect,D] = eig(0.5.*((exp(1i*theta).*A)+(exp(1i*theta).*A)'));
    [~,b] = max(D*ones(m,1));
    V = vect(:,b)./norm(vect(:,b));
    pto(ind) = V'*A*V;
end

output = pto + (tr*ones(1,k)) + (EPS*ones(1,k));
output = [output output(1)];