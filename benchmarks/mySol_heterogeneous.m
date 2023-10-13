function [solU, solDx, solDy, solF, solVx, solVy] = mySol_heterogeneous(x,y,k)

global TAGbench omega;

switch TAGbench
    case 'cavity'
        [solU, solDx, solDy] = cavity(x,y,k);
        solF = 0*x+1;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
    case 'waveguide'
        L = 4.;
        theta = 30.*(pi/180.);
        [solU, solDx, solDy] = waveguide(x,y,k,L,theta);
        solF = 0*x;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
    case 'open'
        theta = pi/4;
        solU  = exp(1i*k*(cos(theta)*x+sin(theta)*y));
        solF  = 0*x;
        solDx = 1i*k*cos(theta) * solU;
        solDy = 1i*k*sin(theta) * solU;
        solVx = cos(theta) * solU;
        solVy = sin(theta) * solU;
    case 'cavity_heterogeneous'
        [solU, solDx, solDy] = cavity(x,y,k);
        solF = 0*x+1;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
    case 'waveguide_heterogeneous'
        L = 4.;
        theta = 30.*(pi/180.);
        [solU, solDx, solDy] = waveguide(x,y,k,L,theta);
        solF = 0*x;
        solVx = solDx/(1i*k);
        solVy = solDy/(1i*k);
    case 'open_heterogeneous'
        k = omega / 1;
        theta = pi/4;
        solU  = exp(1i*k*(cos(theta)*x+sin(theta)*y));
        solF  = 0*x;
        solDx = 1i*k*cos(theta) * solU;
        solDy = 1i*k*sin(theta) * solU;
        solVx = cos(theta) * solU;
        solVy = sin(theta) * solU;
    otherwise
        warning('Error - No valid benchmark has been set.')
end

end