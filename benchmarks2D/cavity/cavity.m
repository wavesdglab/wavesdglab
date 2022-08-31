function [u,dx,dy] = cavity(x,y,k)

N  = 100;

u  = zeros(size(x));
dx = zeros(size(x));
dy = zeros(size(x));
for i=1:2:N
    sinI = sin(i*pi*x);
    cosI = cos(i*pi*x);
    for j=1:2:N
        sinJ = sin(j*pi*y);
        cosJ = cos(j*pi*y);
        c = 16./(i*pi*i*pi + j*pi*j*pi - k*k)/(i*pi*j*pi);
        u  = u  +        c * sinI .* sinJ;
        dx = dx + i*pi * c * cosI .* sinJ;
        dy = dy + j*pi * c * sinI .* cosJ;
    end
end

end