function displayErrors(errors, M, err_treshold, priority_treshold, save_plots, plot_name, method)
    svalues = [];
    evalues = [];
    for k = 1:length(M.s)-1
        svalues = [svalues, M.s(k), M.s(k+1)];
        evalues = [evalues, errors(k), errors(k)];
    end

    f = figure(5); clf; if ~method.skip_plot_position; f.Position = [15.4000   82.6000  500.0000  300.0000]; end
    semilogy(svalues, evalues, '-b', 'linewidth', 1); hold on;
    semilogy((M.s(2:end)+M.s(1:end-1))./2, errors, '.b', 'linewidth', 1);
    xlabel('$s$', 'interpreter', 'latex');
    ylabel('$\epsilon_{rel}$', 'interpreter', 'latex');
    xline(M.s, ':k');
    xlim([M.s(1), M.s(end)]);
    ylim([0.1*min(errors), 10*max([errors; err_treshold])]);
    yline(err_treshold, 'r', 'linewidth', 1);
    %yline(priority_treshold, 'color', [1, 0.7, 0.7], 'linewidth', 1);
    
    if save_plots
        saveas(gca, plot_name+"_errors.eps", 'epsc');
        saveas(gca, plot_name+"_errors.fig", 'fig');
        saveas(gca, plot_name+"_errors.png", 'png');
    end
end