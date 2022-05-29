function displayQualityMetricsSeparately(quality_metrics, save_plots, plot_name, method)   
    ylabels = {'$t_f$', '$\# variables$', '$log_{10}(\epsilon_{max})$', '$feasability$', '$geometric$ $difference$'};
    metric = ["tf", "nbVars", "err", "feasibility"];
    quality_metrics(3,:) = log10(quality_metrics(3,:));
    for i = 1:min(length(ylabels), 4)
        f = figure(6+i); clf; if ~method.skip_plot_position; f.Position = [1.02e+03,1.034e+02,500,300]; end
        plot(quality_metrics(i,:), '.-b', 'linewidth', 1); hold on;
        xlabel('$iteration$', 'interpreter', 'latex');
        ylabel(ylabels{i}, 'interpreter', 'latex');
        xticks(1:size(quality_metrics,2));
        grid on;
        if save_plots
            saveas(gca, plot_name+"_quality_metrics_"+metric(i)+".eps", 'epsc');
            saveas(gca, plot_name+"_quality_metrics_"+metric(i)+".fig", 'fig');
            saveas(gca, plot_name+"_quality_metrics_"+metric(i)+".png", 'png');
        end
    end
end