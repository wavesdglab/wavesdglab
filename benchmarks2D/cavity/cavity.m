function [u,dx,dy] = cavity(x,y,k)

N   = 150;
tol = 1e-8;

num = 1:2:N;
sinI = zeros(size(x(:),1),size(num,2));
cosI = zeros(size(x(:),1),size(num,2));
sinJ = zeros(size(y(:),1),size(num,2));
cosJ = zeros(size(y(:),1),size(num,2));
for n=1:size(num,2)
    i = num(n);
    sinI(:,n) = sin(i*pi*x);
    cosI(:,n) = cos(i*pi*x);
    sinJ(:,n) = sin(i*pi*y);
    cosJ(:,n) = cos(i*pi*y);
end

u  = zeros(size(x));
dx = zeros(size(x));
dy = zeros(size(x));
for m=1:size(num,2)
    for n=1:size(num,2)
        i = num(m);
        j = num(n);
        c = 16./(i*pi*i*pi + j*pi*j*pi - k*k)/(i*pi*j*pi);
        if(abs(c) > tol)
            u(:)  = u(:)  +        c * sinI(:,m) .* sinJ(:,n);
            dx(:) = dx(:) + i*pi * c * cosI(:,m) .* sinJ(:,n);
            dy(:) = dy(:) + j*pi * c * sinI(:,m) .* cosJ(:,n);
        end
    end
end

% u  = zeros(size(x));
% dx = zeros(size(x));
% dy = zeros(size(x));
% for i=1:2:N
%     sinI = sin(i*pi*x);
%     cosI = cos(i*pi*x);
%     for j=1:2:N
%         c = 16./(i*pi*i*pi + j*pi*j*pi - k*k)/(i*pi*j*pi);
%         if(abs(c) > tol)
%             sinJ = sin(j*pi*y);
%             cosJ = cos(j*pi*y);
%             u  = u  +        c * sinI .* sinJ;
%             dx = dx + i*pi * c * cosI .* sinJ;
%             dy = dy + j*pi * c * sinI .* cosJ;
%         end
%     end
% end


end

% N   = 150;
% tol = 1e-9;
% k   = 10*pi;
% 
% MAT = zeros(N,N);
% for i=1:2:N
%     for j=1:2:N
%         MAT(i,j) = 16./(i*pi*i*pi + j*pi*j*pi - k*k)/(i*pi*j*pi);
%     end
% end
% imagesc(MAT)

% IDA = 0;
% IDB = 0;
% for i=1:2:N
%     for j=1:2:N
%         IDA = IDA+1;
%         c = 16./(i*pi*i*pi + j*pi*j*pi - k*k)/(i*pi*j*pi);
%         if(abs(c) > tol)
%             IDB = IDB+1;
%         end
%     end
% end
% [IDA IDB]