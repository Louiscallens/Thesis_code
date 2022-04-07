function displayQualityMetrics(quality_metrics, method)
    figure(4); clf;
    
    ylabels = {'$t_f$', '$geometric$ $difference$', '$log_{10}(\epsilon_{max})$', '$feasability$'};
    for i = 1:4
        subplot(2,2,i);
        plot(quality_metrics(i,:), '.-b');
        xlabel('$iteration$', 'interpreter', 'latex');
        ylabel(ylabels{i}, 'interpreter', 'latex');
        curtick = get(gca, 'xTick');
        xticks(unique(round(curtick)));
    end
    
    if method.save_plots
        saveas(gca, method.plot_name+"_quality_metrics.eps", 'epsc');
        saveas(gca, method.plot_name+"_quality_metrics.fig", 'fig');
        saveas(gca, method.plot_name+"_quality_metrics.png", 'png');
    end
end