function [solU, solDx, solDy, solF, solVx, solVy] = mySol(x,y)

global k TAGbench;

switch TAGbench
    case 'cavity'
        [solU, solDx, solDy] = cavity(x,y,k);
        solF = 0*x+1;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
%         A = abs(x(1)-x(end));
%         if A~=0
%             A = A/abs(A);
%         end
%         B = abs(y(1)-y(end));
%         if B~=0
%             B = B/abs(B);
%         end
%         solDsDs = A * solDx + B * solDy;
    case 'waveguide'
        L = 1.;
        L = 4.;
        theta = 30.*(pi/180.);
        [solU, solDx, solDy] = waveguide(x,y,k,L,theta);
        solF = 0*x;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
    case 'open'
        theta = pi/4; %pi/4
        solU  = exp(1i*k*(cos(theta)*x+sin(theta)*y));
        solF  = 0*x;
        solDx = 1i*k*cos(theta) * solU;
        solDy = 1i*k*sin(theta) * solU;
        solVx = cos(theta) * solU;
        solVy = sin(theta) * solU;
%         A = sign(abs(x(1)-x(end)));
%         B = sign(abs(y(1)-y(end)));
%         solDsDs = - (k^2) * (cos(theta) * A + sin(theta) * B)^2 * solU;
    otherwise
        warning('Error - No valid benchmark has been set.')
end
end