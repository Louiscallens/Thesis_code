function I = get_integration_matrix(tau)
    N = length(tau);
    I = zeros(N);
    for i = 1:N
        for j = 1:N
            I(i,j) = integral(@(t)evaluateLagrangian(tau, j, t), tau(1), tau(i));
        end
    end
end

function res = evaluateLagrangian(tau, i, t)
    Y = 0.*tau;
    Y(i) = 1;
    res = LagrangePolynomialEval(tau, Y, t);
end