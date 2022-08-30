function [u,dx,dy] = cavity(x,y,k)
    N=100;
    
    u   = zeros(size(x));
    dx  = zeros(size(x));
    dy  = zeros(size(x));
    
%     for i=1:2:N
%         for j=1:2:N
%             ipi = i*pi;
%             jpi = j*pi;
%             
%             ipi2 = ipi*ipi;
%             jpi2 = jpi*jpi;
%             
%             c = 16./(ipi2+jpi2-k*k)/(ipi*jpi);
%             
%             u   = u    +      c*sin(ipi*x).*sin(jpi*y);
%             
%             dx  = dx   + ipi *c*cos(ipi*x).*sin(jpi*y);
%             dy  = dy   + jpi *c*sin(ipi*x).*cos(jpi*y);
%         end
%     end

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
