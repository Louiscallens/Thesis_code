function res = LagrangePolynomialDiff(X,Y,x) 
    x_sym = casadi.MX.sym('x');
    f = casadi.Function('f', {x_sym}, {LagrangePolynomialEval(X, Y, x_sym)});
    g = casadi.Function('g', {x_sym}, {jacobian(f(x_sym), x_sym)});
    res = full(g(x));
end


