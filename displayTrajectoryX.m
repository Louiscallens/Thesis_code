function displayTrajectoryX(res, M, problem, save_plots, plot_name, method)
    svalues = [];
    tvalues = [];
    tfine = [];
    for k = 1:length(res.tc)
        tfine = [tfine, linspace(res.tc{k}(1), res.tc{k}(end), 100)];
        tvalues = [tvalues, res.tc{k}(1:end-1)]; 
        svalues = [svalues, M.sc{k}(1:end-1)];
    end
    tvalues = [tvalues, res.t(end)];
    svalues = [svalues, M.s(end)];
    
    xfine = [];
    xvalues = [];
    xvaluesfirst = [];
    for k = 1:length(res.X)
        if k < length(res.X); xfine = [xfine, LagrangePolynomialEval(res.tc{k}, [res.X{k}, res.X{k+1}(:,1)], linspace(res.tc{k}(1), res.tc{k}(end), 100))]; end
        xvalues = [xvalues, res.X{k}]; 
        xvaluesfirst = [xvaluesfirst, res.X{k}(:,1)];
    end
    
    [tvalues_u, uvalues, tVars, uVars] = get_uvalues(res, M);
    if size(uVars,1) == 1
        uVars(2,:) = zeros(size(uVars(1,:)));
        uvalues(2,:) = zeros(size(uvalues(1,:)));
    end
    
    ufine = interp1(tvalues_u, uvalues(2,:), tfine);
    
    s_marks = 100:100:M.s(end);
    t_marks = NaN + zeros(size(s_marks));
    for i = 1:length(s_marks)
        t_marks(i) = integral(@(time) 1./get_s_derivative(problem.myTrack, interp1(svalues, [res.X{:}]', time)', time), 0, s_marks(i), 'arrayvalued', 1);
    end
    
    f = figure(1); clf; if ~method.skip_plot_position; f.Position = [25.0000  464.2000  500.0000  300.0000]; end
    ylabels = {'$e(t)$', '$\psi(t)$', 'v(t)'};
    
    for i = 1:problem.nx
        subplot(problem.nx, 1, i);
        try xline(t_marks, '-', 'color', [0.8, 0.8, 0.8], 'linewidth', 1); catch; end; hold on;
        plot(tvalues, xvalues(i,:), '.-', 'linewidth', 1); hold on;
        %plot(tfine, xfine(i,:), '-b', 'linewidth', 1); hold on; plot(tvalues, xvalues(i,:), '.b', 'linewidth', 1);
        plot(res.t, xvaluesfirst(i,:), '.k', 'linewidth', 1);
        %xline(res.t, ':', 'color', [0.5, 0.5, 0.5]);
        xlabel('$t$', 'interpreter', 'latex');
        ylabel(ylabels{i}, 'interpreter', 'latex');
        xlim([res.t(1), res.t(end)]);
        if i == 1
            plot(tvalues,  problem.b(tvalues), '-r', 'linewidth', 1);
            plot(tvalues, -problem.b(tvalues), '-r', 'linewidth', 1);
        elseif i == 3
            plot(tVars, problem.max_v.*problem.roll_off(uVars(2,:)), '.-r', 'linewidth', 1);
            %plot(tfine, problem.max_v.*problem.roll_off(ufine), '-r', 'linewidth', 1); hold on; plot(tVars, problem.max_v.*problem.roll_off(uVars(2,:)), '.r', 'linewidth', 1);
        end
    end
    
    if save_plots
        %saveas(gca, "figs/case_studies/"+plot_name+"_all_states.eps", 'epsc');
        %saveas(gca, "figs/case_studies/"+plot_name+"_all_states.fig", 'fig');
        %saveas(gca, "figs/case_studies/"+plot_name+"_all_states.png", 'png');
        saveas(gca, plot_name+"_all_states.eps", 'epsc');
        saveas(gca, plot_name+"_all_states.fig", 'fig');
        saveas(gca, plot_name+"_all_states.png", 'png');
    end
    %{
    xlim([0.15, 0.9]);
    %ylim([47.7, 49.7]);
    ylim([47.7, 50.1]);
    if save_plots
        saveas(gca, plot_name+"_all_states_zoom1.eps", 'epsc');
        saveas(gca, plot_name+"_all_states_zoom1.fig", 'fig');
        saveas(gca, plot_name+"_all_states_zoom1.png", 'png');
    end
    
    xlim([2.4, 3.8]); 
    %ylim([48.5706, 48.5716]);
    ylim([48.57115, 48.57135]);
    if save_plots
        saveas(gca, plot_name+"_all_states_zoom2.eps", 'epsc');
        saveas(gca, plot_name+"_all_states_zoom2.fig", 'fig');
        saveas(gca, plot_name+"_all_states_zoom2.png", 'png');
    end
    %}
    
    y = [res.X{:}];
    f = figure(2); clf; if ~method.skip_plot_position; f.Position = [1041,4.698e+02,500,3e+02]; end
    fine_svalues = linspace(0, M.s(end), 1000);
    fine_angles = problem.myTrack.evaluate_angle(fine_svalues);
    angles = problem.myTrack.evaluate_angle(svalues);
    [fine_Xs, fine_Ys] = problem.myTrack.evaluate_track_param(fine_svalues);
    [Xs, Ys] = problem.myTrack.evaluate_track_param(svalues);
    plot(fine_Xs, fine_Ys, '--k'); hold on; %center-line
    plot(fine_Xs-problem.b(fine_svalues).*sin(fine_angles), fine_Ys+problem.b(fine_svalues).*cos(fine_angles), '-k', 'linewidth', 1); % left bound
    plot(fine_Xs+problem.b(fine_svalues).*sin(fine_angles), fine_Ys-problem.b(fine_svalues).*cos(fine_angles), '-k', 'linewidth', 1); % right bound
    plot(Xs - y(1,:).*sin(angles), Ys + y(1,:).*cos(angles), '.-r', 'linewidth', 1);
    
    y1 = xvaluesfirst;
    angles1 = problem.myTrack.evaluate_angle(M.s);
    [Xs1, Ys1] = problem.myTrack.evaluate_track_param(M.s);
    plot(Xs1 - y1(1,:).*sin(angles1), Ys1 + y1(1,:).*cos(angles1), '.k', 'linewidth', 1);
    axis equal;
    
    angle_marks = problem.myTrack.evaluate_angle(s_marks);
    [Xs_marks, Ys_marks] = problem.myTrack.evaluate_track_param(s_marks);
    for i = 1:length(s_marks)
        plot([Xs_marks(i)+problem.b(s_marks(i))*sin(angle_marks(i)), Xs_marks(i)-problem.b(s_marks(i))*sin(angle_marks(i))], ...
             [Ys_marks(i)-problem.b(s_marks(i))*cos(angle_marks(i)), Ys_marks(i)+problem.b(s_marks(i))*cos(angle_marks(i))], 'color', [0.8, 0.8, 0.8], 'linewidth', 1);
    end
    plot([0,0], [-problem.b(0), problem.b(0)], 'k', 'linewidth', 1);
    
    if save_plots
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajX.eps", 'epsc');
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajX.fig", 'fig');
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajX.png", 'png');
        saveas(gca, plot_name+"_trajX.eps", 'epsc');
        saveas(gca, plot_name+"_trajX.fig", 'fig');
        saveas(gca, plot_name+"_trajX.png", 'png');
    end
    
    %xlim([8, 17]);
    %ylim([3.5, 10]);
    
    %if save_plots
    %    saveas(gca, plot_name+"_trajX_zoomed.eps", 'epsc');
    %    saveas(gca, plot_name+"_trajX_zoomed.fig", 'fig');
    %    saveas(gca, plot_name+"_trajX_zoomed.png", 'png');
    %end
    
end