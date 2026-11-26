% Assumption: x in [-1,1]

function val = functionsShapeDerLIN(x,degree)

global Options
switch Options.Basis

    case 'Bernstein'
        val = functionsBernsteinDer(x,degree);
    case 'Lagrange'
        val = functionsLagrangeDerLIN(x,degree);
    case 'Lobbato'
        val = functionsLobbatoDer(x,degree);

    otherwise
        N = degree+1;
        x = x(:);
        val = zeros(size(x,1),N);
        functionsKernel = functionsJacobi(x,1,1,N-2);
        functionsKernelDer = functionsJacobiDer(x,1,1,N-2);

        % nodal modes
        val(:,1) = -0.5;
        val(:,2) =  0.5;

        % edge modes
        for n=1:N-2
            val(:,n+2) = 0.25 * (1-x) .* (1+x) .* functionsKernelDer(:,n) ...
                - 0.5 * x .* functionsKernel(:,n);
        end
end

end