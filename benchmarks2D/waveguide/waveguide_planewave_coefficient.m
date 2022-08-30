function g = waveguide_planewave_coefficient(n,k,theta)
threshold = 1.e-10;

npi = n*pi;

a = k*sin(theta);

ip = 1i*(a+npi);
im = 1i*(a-npi);

Ip = 1.;
Im = 1.;

if (abs(ip) > threshold)
    Ip = (exp(ip)-1.)/ip;
end

if (abs(im) > threshold)
    Im = (exp(im)-1.)/im;
end

g = (Ip-Im)/(2.*1i);
end