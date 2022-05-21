function displayTrajectoryX(res, M, problem, save_plots, plot_name)
    svalues = [];
    tvalues = [];
    for k = 1:length(res.tc)
        %if k == 1
        %   tvalues = [tvalues, linspace(res.tc{k}(1), res.tc{k+1}(1), 100)];
        %else
           tvalues = [tvalues, res.tc{k}(1:end-1)]; 
        %end
        svalues = [svalues, M.sc{k}(1:end-1)];
    end
    tvalues = [tvalues, res.t(end)];
    svalues = [svalues, M.s(end)];
    
    xvalues = [];
    xvaluesfirst = [];
    for k = 1:length(res.X)
        %if k == 1
        %   xvalues = [xvalues, LagrangePolynomialEval(res.tc{k}, [res.X{k}, res.X{k+1}(:,1)], linspace(res.tc{k}(1), res.tc{k+1}(1), 100))];
        %else
           xvalues = [xvalues, res.X{k}]; 
        %end
        xvaluesfirst = [xvaluesfirst, res.X{k}(:,1)];
    end
    
    [tvalues_u, uvalues, tVars, uVars] = get_uvalues(res, M);
    if size(uVars,1) == 1
        uVars(2,:) = zeros(size(uVars(1,:)));
        uvalues(2,:) = zeros(size(uvalues(1,:)));
    end
    
    f = figure(1); clf; f.Position = [25.0000  464.2000  500.0000  300.0000];
    ylabels = {'$e(t)$', '$\psi(t)$', 'v(t)'};
    for i = 1:problem.nx
        subplot(problem.nx, 1, i);
        plot(tvalues, xvalues(i,:), '.-', 'linewidth', 1); hold on;
        plot(res.t, xvaluesfirst(i,:), '.k', 'linewidth', 1);
        xline(res.t, ':', 'color', [0.5, 0.5, 0.5]);
        xlabel('$t$', 'interpreter', 'latex');
        ylabel(ylabels{i}, 'interpreter', 'latex');
        xlim([res.t(1), res.t(end)]);
        if i == 1
            plot(tvalues,  problem.b(tvalues), '-r', 'linewidth', 1);
            plot(tvalues, -problem.b(tvalues), '-r', 'linewidth', 1);
        elseif i == 3
            plot(tVars, problem.max_v.*problem.roll_off(uVars(2,:)), '.-r', 'linewidth', 1);
            %plot(tvalues_u, problem.max_v.*problem.roll_off(uvalues(2,:)), '-r', 'linewidth', 1);
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
    
    y = [res.X{:}];
    f = figure(2); clf; f.Position = [1041,4.698e+02,500,3e+02];
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