% Copyright (C) 2023, CNRS, Inria, ENSTA Paris
% See the LICENSE.txt file in the root directory for license information
% Author: Axel Modave

function normal = getNormalTRI(V1,V2,V3)

T1 = V2-V1;
T2 = V3-V2;
T3 = V1-V3;
N1 = T3 - (T3*T1') * T1/norm(T1)^2;
N2 = T1 - (T1*T2') * T2/norm(T2)^2;
N3 = T2 - (T2*T3') * T3/norm(T3)^2;
N1 = N1/norm(N1);
N2 = N2/norm(N2);
N3 = N3/norm(N3);
normal = [N1 ; N2 ; N3];

% hold off
% plot([V1(1) V2(1) V3(1) V1(1)], [V1(2) V2(2) V3(2) V1(2)]);
% hold on
% plot([(V1(1)+V2(1))/2 (V1(1)+V2(1))/2+N1(1)], [(V1(2)+V2(2))/2 (V1(2)+V2(2))/2+N1(2)]);
% plot([(V2(1)+V3(1))/2 (V2(1)+V3(1))/2+N2(1)], [(V2(2)+V3(2))/2 (V2(2)+V3(2))/2+N2(2)]);
% plot([(V3(1)+V1(1))/2 (V3(1)+V1(1))/2+N3(1)], [(V3(2)+V1(2))/2 (V3(2)+V1(2))/2+N3(2)]);

end