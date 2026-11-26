% Assumption: x in [-1,1]

function val = functionsShapeLIN(x,degree)

global Options
switch Options.Basis

    case 'Bernstein'
        val = functionsBernstein(x,degree);
    case 'Lagrange'
        val = functionsLagrangeLIN(x,degree);
    case 'Lobbato'
        val = functionsLobbato(x,degree);

    otherwise % Jacobi
        N = degree+1;
        x = x(:);
        val = zeros(size(x,1),N);
        functionsKernel = functionsJacobi(x,1,1,N-2);

        % nodal modes
        val(:,1) = (1-x)/2;
        val(:,2) = (1+x)/2;

        % edge modes
        for n=1:N-2
            val(:,n+2) = val(:,1) .* val(:,2) .* functionsKernel(:,n);
        end
end

end