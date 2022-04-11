import casadi.*

x0 = 0;

nx = 1;
t = casadi.MX.sym('t');

f = @(t, x) pi.*cos(pi.*t);
x_e = @(t) x0 + sin(pi.*t);
%f = @(t, x) 0.*(t < -1/2) + pi.*cos(pi.*t).*(t >= -1/2).*(t < 1/2);
%x_e = @(t) x0.*(t < -1/2) + (x0 + 1 + sin(pi.*t)).*(t >= -1/2).*(t < 1/2) + (x0+2).*(t>1/2);

err_ests = [];
err_real = [];

for Nk = 4:4:24
    X0 = casadi.MX.sym('X0');
    Xcoll = casadi.MX.sym('Xcoll', Nk, 1);
    tau = get_collocation_interval(-1, 1, Nk);
    polEval = casadi.Function('polEval', {X0, Xcoll, t}, {LagrangePolynomialEval(tau, [X0; Xcoll]', t)});
    dot_Pi = casadi.Function('dot_Pi', {X0, Xcoll, t}, {jacobian(polEval(X0, Xcoll, t), t)});

    g = {};
    g{1} = dot_Pi(X0, Xcoll, tau(1)) - f(tau(1), X0);
    for j = 2:Nk
        g{end+1} = dot_Pi(X0, Xcoll, tau(j)) - f(tau(j), Xcoll(j-1));
    end
    g = vertcat(g{:});
    
    rf = rootfinder('rf', 'newton', struct('x', Xcoll, 'p', X0, 'g', g));
    res = rf('x0', repmat(x0, Nk, 1), 'p', x0);
    results = full(res.x);
    results = [x0; results];
    
    
    [est, ~] = get_coll_error_est(tau, results', results(1:end-1)', 0, @(x, u, t)f(t, x));
    err_ests = [err_ests, est];
    err_real = [err_real, get_actual_error(tau, results', x_e)];
    
    %pause(1);
    %plot(tvalues, LagrangePolynomialEval(tau, results', tvalues), 'displayname', "Nk = " + num2str(Nk));
    %tvalues = linspace(-1,1,1000);
    %xdvalues = NaN + zeros(size(tvalues));
    %for i = 1:length(tvalues)
    %    xdvalues(i) = full(dot_Pi(results(1), results(2:end), tvalues(i)));
    %end
    %plot(tvalues, xdvalues);
end
%legend();

figure;
scatter(4:4:24, log10(err_ests)); hold on;
scatter(4:4:24, log10(err_real));
legend('est', 'act', 'autoupdate', 'off');
xline([8, 12, 16, 20], ':');
yline([-5, -10], ':');
%yline([-3, -2.5, -2, -1.5, -1, -0.5], ':');



function inf_norm = get_actual_error(tcoll, xcoll, x_e)
    tvalues = linspace(-1,1,10000);
    sol = @(t) LagrangePolynomialEval(tcoll, xcoll, t);
    inf_norm = 0;
    for i = tvalues
        if abs(sol(i) - x_e(i)) > inf_norm
            inf_norm = abs(sol(i) - x_e(i));
        end
    end
end
