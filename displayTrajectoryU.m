function displayTrajectoryU(res, M, problem, save_plots, plot_name)
    plot_polynomials = false;
    
    tvalues = [];
    uvalues = [];
    tVars = 0;
    uVars = res.U{1};
    for k = 1:length(res.U)-1
        if M.Nu(k) == 0
            tvalues = [tvalues, res.t(k), res.t(k+1)];
            uvalues = [uvalues, res.U{k}, res.U{k}];
            tVars = [tVars, (res.t(k)+res.t(k+1))/2];
            uVars = [uVars, res.U{k}];
        elseif M.Nu(k) == 1
            uvalues = [uvalues, res.U{k} + (res.tc{k}(1:end-1)-res.tc{k}(1)).*(res.U{k+1}-res.U{k})./(res.tc{k}(end)-res.tc{k}(1))];
            tvalues = [tvalues, res.tc{k}(1:end-1)];
            tVars = [tVars, res.t(k), res.t(k+1)];
            uVars = [uVars, res.U{k}, res.U{k+1}(:,1)];
        else
            if plot_polynomials
                temp = linspace(res.tc{k}(1), res.tc{k}(end), 100);
                uvalues = [uvalues, LagrangePolynomialEval(res.tc{k}(1:end-1), res.U{k}, temp)];
                tvalues = [tvalues, temp];
                tVars = [tVars, res.tc{k}(1:end-1)];
                uVars = [uVars, res.U{k}];
            else
                uvalues = [uvalues, res.U{k}];
                tvalues = [tvalues, res.tc{k}(1:end-1)];
                tVars = [tVars, res.tc{k}(1:end-1)];
                uVars = [uVars, res.U{k}];
            end
        end
    end
    uvalues = [uvalues, res.U{end}];
    tvalues = [tvalues, res.t(end)];
    
    max_accel = 20;
    min_accel = -5;
    
    figure(3); clf;
    subplot(211);
    plot(tVars, max_accel.*(1.1.*cos(uVars(2,:))).^4, '.-r'); hold on;
    plot(tVars, min_accel.*(1.1.*cos(uVars(2,:)).^4), '.-r');
    plot(tvalues, uvalues(1,:), 'b', 'linewidth', 1); hold on;
    plot(tVars, uVars(1,:), '.b');
    xline(res.t, ':', 'color', [0.5, 0.5, 0.5]);
    xlabel('$t$', 'interpreter', 'latex');
    title('$u_1(t)$ (throttle)', 'interpreter', 'latex');
    subplot(212);
    plot(tvalues, uvalues(2,:), 'b', 'linewidth', 1); hold on;
    plot(tVars, uVars(2,:), '.b');
    xline(res.t, ':', 'color', [0.5, 0.5, 0.5]);
    xlabel('$t$', 'interpreter', 'latex');
    title('$u_2(t)$ (steering angle)', 'interpreter', 'latex');
    
    if save_plots
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajU.eps", 'epsc');
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajU.fig", 'fig');
        %saveas(gca, "figs/case_studies/"+plot_name+"_trajU.png", 'png');
        saveas(gca, plot_name+"_trajU.eps", 'epsc');
        saveas(gca, plot_name+"_trajU.fig", 'fig');
        saveas(gca, plot_name+"_trajU.png", 'png');
    end
end