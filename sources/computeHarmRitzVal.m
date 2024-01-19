% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Timothee Raynaud, Axel Modave, Pierre Marchand

% This function computes the HRV (reps. RV) at step i and the distance between the iEig smallest eigenvalues and the HRV (resp. RV)
function [hrv, dist_hrv, rv, dist_rv] = computeHarmRitzVal(H, cs, sn, i, A, iEig)

% HRV and RV
matQ = eye(i+1,i+1);
for j = 1:i
    tempQ = eye(i+1,i+1);
    tempQ(j,j) = cs(j)';
    tempQ(j+1,j+1) = cs(j);
    tempQ(j,j+1) = sn(j)';
    tempQ(j+1,j) = -sn(j)';
    matQ = tempQ*matQ;
end

R = zeros(i+1, i);

R(1:i,1:i) = H(1:i,1:i);

Ht_i = matQ'*R;

H_i = Ht_i(1:i,1:i);

[~,D] = eig(Ht_i'*Ht_i,H_i');
hrv = diag(D);

rv = eig(H_i);

% Distances
dist_hrv = zeros(1,iEig);
dist_rv = zeros(1,iEig);

lambdaMin = eigs(A, iEig, 'smallestabs');

for j = 1:iEig
    dist_hrv(j) = min(abs(hrv - lambdaMin(j)));
    dist_rv(j) = min(abs(rv - lambdaMin(j)));
end

end
