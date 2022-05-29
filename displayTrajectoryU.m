function displayTrajectoryU(res, M, problem, save_plots, plot_name, method)
       
    [tvalues, uvalues, tVars, uVars] = get_uvalues(res, M);
    if size(uVars,1) == 1
        uVars(2,:) = zeros(size(uVars(1,:)));
        uvalues(2,:) = zeros(size(uvalues(1,:)));
    end
    
    f = figure(3); clf; if ~method.skip_plot_position; f.Position = [536.2000  467.4000  500.0000  300.0000]; end
    subplot(211);
    plot(tVars, problem.max_accel.*problem.roll_off(uVars(2,:)), '.-r'); hold on;
    plot(tVars, problem.min_accel.*problem.roll_off(uVars(2,:)), '.-r');
    %plot(tvalues, problem.max_accel.*problem.roll_off(uvalues(2,:)), '-r', 'linewidth', 1);
    %plot(tvalues, problem.min_accel.*problem.roll_off(uvalues(2,:)), '-r', 'linewidth', 1);
    %plot(tvalues, uvalues(1,:), 'b', 'linewidth', 1); hold on;
    plot(tVars, uVars(1,:), '.-b', 'linewidth', 1);
    xline(res.t, ':', 'color', [0.5, 0.5, 0.5]);
    xlim([res.t(1), res.t(end)]);
    yline(0, 'color', [0.7, 0.7, 0.7]);
    ylim([problem.min_accel-2, problem.max_accel+2]);
    xlabel('$t$', 'interpreter', 'latex');
    ylabel('$u_1(t)$', 'interpreter', 'latex');
    
    subplot(212);
    %plot(tvalues, uvalues(2,:), 'b', 'linewidth', 1); hold on;
    plot(tVars, uVars(2,:), '.-b', 'linewidth', 1);
    xline(res.t, ':', 'color', [0.5, 0.5, 0.5]);
    xlim([res.t(1), res.t(end)]);
    yline(0, 'color', [0.7, 0.7, 0.7]);
    xlabel('$t$', 'interpreter', 'latex');
    ylabel('$u_2(t)$', 'interpreter', 'latex');
    %}
    if save_plots
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajU.eps", 'epsc');
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajU.fig", 'fig');
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajU.png", 'png');
        saveas(gca, plot_name+"_trajU.eps", 'epsc');
        saveas(gca, plot_name+"_trajU.fig", 'fig');
        saveas(gca, plot_name+"_trajU.png", 'png');
    end
end