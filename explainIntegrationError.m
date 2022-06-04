%main;

%k = 7;
%k = 23;
k = 4;
%k = 40;

track = problem.myTrack;
x = [res.X{k}, res.X{k+1}(:,1)]; u = [res.U{k}, res.U{k+1}(:,1)]; sc = M.sc{7};
sdot = get_s_derivative(track, x, sc);

svals = linspace(sc(1), sc(end), 1000);
figure(11); clf;

rho = @(s) track.evaluate_radius_curvature(s);
x1 = @(s) LagrangePolynomialEval(sc, x(1,:), s); 
x2 = @(s) LagrangePolynomialEval(sc, x(2,:), s);
x3 = @(s) LagrangePolynomialEval(sc, x(3,:), s);
%x1 = @(s) interp1(sc, x(1,:), s); x2 = @(s) interp1(sc, x(2,:), s); x3 = @(s) interp1(sc, x(3,:), s);

% actual intergrand (solution substituted)
g = @(s) 1./(rho(s)./(rho(s)-x1(s)).*x3(s).*cos(x2(s)));

figure(11); clf;
plot(svals, g(svals), '-b', 'linewidth', 1); hold on;
xlim([svals(1), svals(end)]);

% trapezoid approximation
h = @(s) 1./interp1(sc, sdot, s);

plot(svals, h(svals), 'g', 'linewidth', 1);

if ~method.save_plots
    legend('$1/\dot{s}$', 'trapezoid approx', 'location', 'best', 'autoupdate', 'off');
end
scatter(sc, g(sc), 30, 'b', 'filled');

if method.save_plots
    %saveas(gca, method.plot_name+"_integrand.eps", 'epsc');
    %saveas(gca, method.plot_name+"_integrand.fig", 'fig');
    %saveas(gca, method.plot_name+"_integrand.png", 'png');
end

t1 = integral(@(s) g(s), svals(1), svals(end));
t2 = integral(@(s) h(s), svals(1), svals(end));
I = get_integration_matrix(sc);
t3 = transpose(I(end,:)*(1./sdot'));
disp("real time:      "+num2str(t1));
disp("trapezoid time: "+num2str(t2));
disp("coll obj time:  "+num2str(t3));
disp("diff:           "+num2str((t3-t2)));
%figure(2); clf;
%semilogy(svals, f(svals)-g(svals));

