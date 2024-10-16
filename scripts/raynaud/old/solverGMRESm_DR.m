% G. Meurant and J. Duintjer Tebbens, Krylov Methods for Nonsymmetric Linear Systems: From Theory to Computations. Springer International Publishing, 2020. doi: 10.1007/978-3-030-55251-0.


function [x,ni,nc,resn] = solverGMRESm_DR(A,b,x0,epsi,nitmax,m,kR)
%
m = min(m,nitmax);
n = size(A,1);
rhs = zeros(m+1,1);
resn = zeros(1,nitmax+1);
nb = norm(b);
x = x0;
r = b - A * x;
rzero = r;
H = zeros(m+1,m);
HH = H;
bet = norm(r);
nc = 1;% number of cycles
ni = 0;% number of iterations
iconv = 0;
resn(1) = bet;
rhs(1) = bet;
rhsR = rhs;
%
% do a first cycle of GMRES(m)
V = zeros(n,m+1);
rot = zeros(2,m);
v = r / bet;
V(:,1) = v;
for k = 1:m
    ni = ni + 1;% number of iterations
    [V,H,HH,rhs,rhsR,rot] = gmres_iter(k,A,V,H,HH,rhs,rhsR,rot, rzero,0,[]);
    nresidu = abs(rhs(k+1));
    resn(ni+1) = nresidu;
    % convergence test or too many iterations
    if nresidu < (epsi * nb) || ni >= nitmax % convergence 
        iconv = 1;
    break
    end % if nresidu
end % for k - end of the first cycle
% computation of the solution at the end of the first cycle
y = triu(H(1:k,1:k)) \ rhs(1:k);
x = x0 + V(:,1:k) * y;
if iconv == 1
    % we have to stop
    resn = resn(1:ni+1);
    return
end % if iconv
% residual vector of the least squares problem
rR = rhsR - HH(1:k+1,1:k) * y;
% start the while loop with computation of the harmonic Ritz
% values and vectors sorted by modulus
kR0 = kR;
[g,kR] = harm_Ritz_vec(kR0,HH);
QQ = [];
%
while nresidu > (epsi * nb) && (ni < nitmax)
    % —- Loop on cycles 
    nc = nc + 1;% number of cycles
    if kR >= m 
        fprintf('\n Error: kR = %d is too large compare to m = %d, ... stop \n\n',kR,m)
        return
    end
    rot = zeros(2,m);% init Givens rotations
    rhs = zeros(m+1,1);
    rhsR = rhs;
    r = b - A * x;
    rzero = r;
    x0 = x;
    if kR ~=0
        % first phase
        % orthogonalize the harmonic Ritz vectors
        P = orth_mgs(g,'dreorth');
        P = [P;zeros(1,size(P,2))];% append a zero row
        PP = zeros(m+1,kR+1);
        PP(:,1:kR) = P;
        % double orthogonalization of rR against P
        w = orth_vec(P,rR);
        PP(:,kR+1) = w;
        V(:,1:kR+1) = V(:,1:m+1) * PP;% project V 
        % orthogonalize v_(kR+1) against the previous basis vectors 
        w = orth_vec(V(:,1:kR),V(:,kR+1));
        V(:,kR+1) = w;
        H(1:kR+1,1:kR) = PP' * HH(1:m+1,1:m)* P(1:m,:);% project H 
        HH(1:kR+1,1:kR) = H(1:kR+1,1:kR);% save it
        HH(kR+2:m+1,1:kR) = zeros(m-kR,kR);
        H(kR+2:m+1,1:kR) = zeros(m-kR,kR);
        [QQ,RR] = qr(HH(1:kR+1,1:kR));% QR factorization 
        H(1:kR+1,1:kR) = RR(1:kR+1,1:kR);
        % part of the right-hand side 
        rhs(1:kR+1) = V(:,1:kR+1)' * rzero;
        rhsR(1:kR+1) = rhs(1:kR+1);% save it
        rhs(1:kR+1)=QQ' * rhs(1:kR+1);% modify the right-hand side
    else % if kR, case without augmentation 
        rhs(1) = norm(rzero);
        V(:,1) = rzero / rhs(1);
    end % if kR
    % second phase: start a cycle of GMRES(m) from what has been
    % already computed
    for k = kR+1:m
        ni = ni + 1;% number of iterations
        [V,H,HH,rhs,rhsR,rot] = gmres_iter(k,A,V,H,HH,rhs,rhsR,rot,rzero,kR,QQ);
        nresidu = abs(rhs(k+1));
        resn(ni+1) = nresidu;
        % convergence test or too many iterations
        if nresidu < (epsi * nb) || ni >= nitmax 
            % convergence
            iconv = 1;
            break % get out of the loop
        end % if nresidu
    end % for k - end of one cycle 
        % computation of the solution at the end of the cycle
    y = triu(H(1:k,1:k)) \ rhs(1:k);
    x = x0 + V(:,1:k) * y;
    if iconv == 1 % we have to stop
        resn = resn(1:ni+1);
        return 
    end % if iconv 
    % we have not converged yet,compute the harmonic eigenvectors
    % and restart 
    if kR0 == 0
        g = [];
    else
        [g,kR] = harm_Ritz_vec(kR0,HH);
    end % if kR0 
    rR = rhsR - HH(1:k+1,1:k) * y;
end % while, loop on cycles % if we get here we have done the max number of cycles,
    % may be without convergence
    resn = resn(1:ni);
end % function 


function w = orth_vec(P,rR)
    % double orthogonalize rR against the columns of P
    kRP = size(P,2);
    w = rR / norm(rR);
    for j = 1:kRP 
        alpha = P(:,j)' * w;
        w = w - alpha * P(:,j);
    end % for 
    for j = 1:kRP
        alpha = P(:,j)' * w;
        w = w - alpha * P(:,j);
    end % for j 
    w = w / norm(w);
end % function 

function V = orth_mgs(A,dreorth)
    % (double) orthogonalisation of the columns of A
    if nargin == 1 
        dreorth = 'nodreorth';
    end 
    [m,n] = size(A);
    V = zeros(m,n);
    v = A(:,1);
    V(:,1) = v / norm(v);
    for k = 2:n
        v = A(:,k);
        for j = 1:k-1 
            alpha = v' * V(:,j);
            v = v - alpha * V(:,j);
        end % for j 
        if strcmpi(dreorth,'dreorth') == 1
            for j = 1:k-1
                alpha = v' * V(:,j);
                v = v - alpha * V(:,j);
            end % for j 
        end % if 
        nv = norm(v);
        if nv <= 1e-15
            fprintf('\n orth_mgs: Breakdown, step %d, val=%g \n\n',k,nv) 
            %return 
        end % if 
        v = v / norm(v);
        V(:,k) = v;
    end % for k 
end % function 

function [g,kR] = harm_Ritz_vec(kR,HH)
    % computation of the harmonic Ritz vectors 
    if kR == 0
        g = [];
    else
        [theta,g] = harm_Ritz_valvec(HH);
        % keep only the first kR ones (with the pairs of complex conjugate eigenvalues)
        if isreal(theta(kR))
            theta = theta(1:kR);
            g = g(:,1:kR);
        else
            if abs(abs(theta(kR)) - abs(theta(kR+1))) <= 1e-14 
                kR = kR + 1;% increase kR to keep a pair
                theta = theta(1:kR);
                g = g(:,1:kR);
            else 
                theta = theta(1:kR);
                g = g(:,1:kR);
            end % if abs
        end % if isreal
        % we need real vectors 
        jj = 1;
        for j=1:kR 
            if isreal(theta(jj))
                jj = jj + 1;
            else 
                g(:,jj+1) = imag(g(:,jj));
                g(:,jj) = real(g(:,jj));
                jj = jj + 2;
            end % if
            if jj > kR 
                break 
            end 
        end % for j
    end % if kR 
end % function 


function [lambda,g] = harm_Ritz_valvec(H)
    %computes the harmonic Ritz values and eigenvectors of the upper Hessenberg matrix H
    n = size(H,1);
    k = n - 1;
    Hk = H(1:k,1:k);
    ek = zeros(k,1);
    ek(k) = 1;
    Hki = Hk' \ ek;
    Hk(:,k) = Hk(:,k) + H(k+1,k)^2 * Hki;
    [X,D] = eig(full(Hk));
    lam = diag(D);
    [lamb,I] = sort(abs(lam));
    lambda = lam(I);
    g = X(:,I);
end % function 

function [V,H,HH,rhs,rhsR,rot] = gmres_iter(k,A,V,H,HH,rhs, rhsR,rot,rzero,kR,QQ)
    % iteration k of GMRES-DR 
    Av = A * V(:,k);% matrix-vector product 
    w = Av;
    for l = 1:k
        vl = V(:,l);
        gl = vl' * w;
        H(l,k) = gl;
        w = w - gl * vl;
    end % for l 
    dk = w;
    gk = norm(dk);
    dk1 = dk / gk;
    % normalization of the new vector
    H(k+1,k) = gk;
    V(:,k+1) = dk1;
    gk1 = gk;
    rhs(k+1) = V(:,k+1)' * rzero;
    % compute the right-hand side 
    rhsR(k+1) = rhs(k+1);
    % save it
    HH(1:k+1,k) = H(1:k+1,k);
    % save the column of H
    if kR ~= 0 
        H(1:kR+1,k) = QQ' * H(1:kR+1,k);
        % apply the Q factor to the last column 
    end % if kR 
    % apply the preceding Givens rotations to the last column just computed 
    for kk = kR+1:k-1 
        g1 = H(kk,k);
        g2 = H(kk+1,k);
        H(kk+1,k) = -rot(2,kk) * g1 + rot(1,kk) * g2;
        H(kk,k) = rot(1,kk) * g1 + conj(rot(2,kk)) * g2;
    end % for kk 
    % compute, store and apply a new rotation to zero the last term in kth column 
    gk = H(k,k);
    if gk == 0 
        rot(1,k) = 0;
        rot(2,k) = 1;
    elseif gk1 == 0 
        rot(1,k) = 1;
        rot(2,k) = 0;
    else 
        cs = sqrt(abs(gk1)^2 + abs(gk)^2);
        if abs(gk) < abs(gk1)
            mu = gk / gk1;
            tau = conj(mu) / abs(mu);
        else 
            mu = gk1 / gk;
            tau = mu / abs(mu);
        end % if 
        % store the rotation for the next columns
        rot(1,k) = abs(gk) / cs;% cosine
        rot(2,k) = abs(gk1) * tau / cs;% sine
    end % if gk % modify the diagonal entry and the right-hand side 
    H(k,k) = rot(1,k) * gk + conj(rot(2,k)) * gk1;
    c = rhs(k);
    rhs(k) = rot(1,k) * c;
    rhs(k+1) = -rot(2,k) * c;
    H(k+1,k) = 0;
end % function 