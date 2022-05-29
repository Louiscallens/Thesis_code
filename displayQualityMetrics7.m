function displayQualityMetrics7(quality_metrics, save_plots, plot_name, method)
    f = figure(4); clf; if ~method.skip_plot_position; f.Position = [1.02e+03,1.034e+02,500,300]; end
    
    ylabels = {'$t_f$', '$\# variables$', '$log_{10}(\epsilon_{max})$', '$feasability$', '$geometric$ $difference$'};
    quality_metrics(3,:) = log10(quality_metrics(3,:));
    i = 1;
    %for i = 1:length(ylabels)
    %    subplot(2,2,i);
        plot(quality_metrics(i,:), '.-b'); hold on              ;
        xlabel('$iteration$', 'interpreter', 'latex');
        ylabel(ylabels{i}, 'interpreter', 'latex');
        xticks(1:size(quality_metrics,2));
        grid on;
        %if i == 1
        %    plot(quality_metrics(end,:), '-r', 'linewidth', 1);
        %end
    %end
    
    if save_plots
        saveas(gca, plot_name+"_quality_metrics.eps", 'epsc');
        saveas(gca, plot_name+"_quality_metrics.fig", 'fig');
        saveas(gca, plot_name+"_quality_metrics.png", 'png');
    end
end