function [u,dx,dy] = waveguide(x,y,k,L,theta)

N = 400;
tol = 1.e-10;

u  = zeros(size(x));
dx = zeros(size(x));
dy = zeros(size(x));

for n=1:N
    
    % COEF
    Ip = 1.;
    Im = 1.;
    ip = k*sin(theta)+n*pi;
    im = k*sin(theta)-n*pi;
    if (abs(ip) > tol)
        Ip = -1i*(exp(1i*ip)-1.)/ip;
    end
    if (abs(im) > tol)
        Im = -1i*(exp(1i*im)-1.)/im;
    end
    g = (Ip-Im)/(2.*1i);
    
    % ODE
    kappa = n*n*pi*pi-k*k;
    v = g * x/(1-1i*k*L);
    d = g * 1/(1-1i*k*L);
    if (abs(kappa) > tol)
        s = sqrt(kappa);
        den1 = (s-1i*k) * exp(s*(L-x)) + (s+1i*k) * exp(-s*(L+x));
        den2 = (s-1i*k) * exp(s*(L+x)) + (s+1i*k) * exp(-s*(L-x));
        v = g   * (1./den1 - 1./den2);
        d = g*s * (1./den1 + 1./den2);
        if isnan(v)
            error('NAN');
        end
    end
    
    % END
    npi = n*pi;
    u  = u  +     v.*sin(npi*y);
    dx = dx +     d.*sin(npi*y);
    dy = dy + npi*v.*cos(npi*y);
end


end