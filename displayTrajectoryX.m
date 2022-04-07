function displayTrajectoryX(res, M, problem, save_plots, plot_name)
    svalues = [];
    tvalues = [];
    for k = 1:length(res.tc)
        tvalues = [tvalues, res.tc{k}(1:end-1)];
        svalues = [svalues, M.sc{k}(1:end-1)];
    end
    tvalues = [tvalues, res.t(end)];
    svalues = [svalues, M.s(end)];
    
    xvalues = [];
    xvaluesfirst = [];
    for k = 1:length(res.X)
        xvalues = [xvalues, res.X{k}];
        xvaluesfirst = [xvaluesfirst, res.X{k}(:,1)];
    end
    
    figure(1); clf;
    for i = 1:problem.nx
        subplot(problem.nx, 1, i);
        plot(tvalues, xvalues(i,:), '.-', 'linewidth', 1); hold on;
        plot(res.t, xvaluesfirst(i,:), '.k', 'linewidth', 1);
        xline(res.t, ':', 'color', [0.5, 0.5, 0.5]);
        xlabel('$t$', 'interpreter', 'latex');
        ylabel("$x_"+num2str(i)+"$", 'interpreter', 'latex');
    end
    
    if save_plots
        %saveas(gca, "figs/case_studies/"+plot_name+"_all_states.eps", 'epsc');
        %saveas(gca, "figs/case_studies/"+plot_name+"_all_states.fig", 'fig');
        %saveas(gca, "figs/case_studies/"+plot_name+"_all_states.png", 'png');
        saveas(gca, plot_name+"_all_states.eps", 'epsc');
        saveas(gca, plot_name+"_all_states.fig", 'fig');
        saveas(gca, plot_name+"_all_states.png", 'png');
    end
    
    y = [res.X{:}];
    figure(2); clf;
    angles = problem.myTrack.evaluate_angle(svalues);
    [Xs, Ys] = problem.myTrack.evaluate_track_param(svalues);
    plot(Xs, Ys, '--k'); hold on; %center-line
    plot(Xs-problem.b.*sin(angles), Ys+problem.b.*cos(angles), '-k', 'linewidth', 1); % left bound
    plot(Xs+problem.b.*sin(angles), Ys-problem.b.*cos(angles), '-k', 'linewidth', 1); % right bound
    plot(Xs - y(1,:).*sin(angles), Ys + y(1,:).*cos(angles), '.-r', 'linewidth', 1);
    
    y1 = xvaluesfirst;
    angles1 = problem.myTrack.evaluate_angle(M.s);
    [Xs1, Ys1] = problem.myTrack.evaluate_track_param(M.s);
    plot(Xs1 - y1(1,:).*sin(angles1), Ys1 + y1(1,:).*cos(angles1), '.k', 'linewidth', 1);
    axis equal;
    
    if save_plots
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajX.eps", 'epsc');
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajX.fig", 'fig');
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajX.png", 'png');
        saveas(gca, plot_name+"_trajX.eps", 'epsc');
        saveas(gca, plot_name+"_trajX.fig", 'fig');
        saveas(gca, plot_name+"_trajX.png", 'png');
    end
    
end