% Assumption: x,y in [-1,1]

function val = functionsShapeTRI(x,y,degree)

assert(size(x,1) == size(y,1));

L = size(x,1);
N = 3 + 3*(degree-1) + (degree-1)*(degree-2)/2;

% Barycentric coordinates
l1 =  (y+1)/2;
l2 = -(x+y)/2;
l3 =  (x+1)/2;

val = zeros(L,N);
n = 1;

% Nodal modes
val(:,n) = l2; n=n+1;
val(:,n) = l3; n=n+1;
val(:,n) = l1; n=n+1;

% Edge modes
kernel1 = functionsJacobi(l3-l2,1,1,degree-1);
kernel2 = functionsJacobi(l1-l3,1,1,degree-1);
kernel3 = functionsJacobi(l2-l1,1,1,degree-1);
for ne = 1:degree-1
    val(:,n) = l2 .* l3 .* kernel1(:,ne); n=n+1;
end
for ne = 1:degree-1
    val(:,n) = l3 .* l1 .* kernel2(:,ne); n=n+1;
end
for ne = 1:degree-1
    val(:,n) = l1 .* l2 .* kernel3(:,ne); n=n+1;
end

% Face modes
for n1 = 1:degree-1
    for n2 = 1:degree-1-n1
        val(:,n) = l1 .* l2 .* l3 .* kernel1(:,n1) .* kernel3(:,n2); n=n+1;
    end
end


% % edge 1
% Nedg = degree+2; % to modify
% xloc = l3-l2;
% coef = 4*l2.*l3 ./ (1. - xloc.^2);
% valEdg = functionsLobbato(xloc, degree);
% size(valEdg)
% for ne = 1:degree-1
%     val(:,n) = coef .* valEdg(:,2+ne); n=n+1;
% end
% 
% % edge 2
% Nedg = degree+2; % to modify
% xloc = l1-l3;
% coef = 4*l3.*l1 ./ (1. - xloc.^2);
% valEdg = functionsLobbato(xloc, degree);
% for ne = 1:degree-1
%     val(:,n) = coef .* valEdg(:,2+ne); n=n+1;
% end
% 
% % edge 3
% Nedg = degree+2; % to modify
% xloc = l2-l1;
% coef = 4*l1.*l2 ./ (1. - xloc.^2);
% valEdg = functionsLobbato(xloc, degree);
% for ne = 1:degree-1
%     val(:,n) = coef .* valEdg(:,2+ne); n=n+1;
% end
% 
% % face
% for n1 = 1:degree-1
%     for n2 = 1:degree-1-n1
%         val(:,n) = l1 .* l2.^n1 .* l3.^n2; n=n+1;
%     end
% end


end