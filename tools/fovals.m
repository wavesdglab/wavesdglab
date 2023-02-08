function output = fovals(A,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Wilmer Henao    wi-henao@uniandes.edu.co
% Department of Mathematics
% Universidad de los Andes
% Colombia
% 
% The Field of values of a matrix is a convex set in the complex plane
% that contains all eigenvalues of the given matrix, this m-file plots 
% boundary points of this set and the eigenvalues of the given matrix.
% 
% fovals(A,k)
% A = The square matrix
% k = The number of steps (500 by default) you can call FoV, without this argument
% 
% Copyright (c) 2003, Wilmer Henao
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% Source: https://fr.mathworks.com/matlabcentral/fileexchange/4679-field-of-values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EPS = 0.0000000000000000000001*1i;
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