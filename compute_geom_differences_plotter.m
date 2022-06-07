set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',14);

load('workspaces/thesis/final_problem3/geom_diffs');

save_plots = false;

% plot all diffs
figure;
semilogy(areas_slack, '.-b', 'linewidth', 1); hold on;
semilogy(areas_hp,    '.-r', 'linewidth', 1);
%plot(areas_ref,   '.-k', 'linewidth', 1);
for i = [1, 2, 3, 5, 8, 10, 15]
    yline(areas_ref(i));
    text(5.2, areas_ref(i)*1.15, "$N = "+num2str(100*i)+"$", 'interpreter', 'latex');
end
xlim([1,length(areas_slack)]);
xlabel('$iterations$', 'interpreter', 'latex');
ylabel('$area [m^2]$', 'interpreter', 'latex');
grid on;
if save_plots
    saveas(gca, "E:/1. Unief/thesis/thesis_latex/figs/geom_diffs/areas.eps", 'epsc'); 
    saveas(gca, "E:/1. Unief/thesis/thesis_latex/figs/geom_diffs/areas.fig", 'fig'); 
end

% plot maximal error
figure;
semilogy(max_diffs_slack, '.-b', 'linewidth', 1); hold on;
semilogy(max_diffs_hp,    '.-r', 'linewidth', 1);
%plot(areas_ref,   '.-k', 'linewidth', 1);
for i = [1, 2, 3, 5, 8, 10, 15]
    yline(max_diffs_ref(i));
    text(5.2, max_diffs_ref(i)*1.15, "$N = "+num2str(100*i)+"$", 'interpreter', 'latex');
end
xlim([1,length(areas_slack)]);
xlabel('$iterations$', 'interpreter', 'latex');
ylabel('$distance [m]$', 'interpreter', 'latex');
grid on;
if save_plots
    saveas(gca, "E:/1. Unief/thesis/thesis_latex/figs/geom_diffs/max_dist.eps", 'epsc'); 
    saveas(gca, "E:/1. Unief/thesis/thesis_latex/figs/geom_diffs/max_dist.fig", 'fig'); 
end


% plot where the errors are
figure;
semilogy(svals_all_slack{end}, fvals_all_slack{end}, 'b', 'linewidth', 1); hold on;
semilogy(svals_all_hp{end}, fvals_all_hp{end}, 'r', 'linewidth', 1);
%plot(svals_all_ref{end}, fvals_all_ref{end}, 'k', 'linewidth', 1);
grid on;
xlabel('$s$', 'interpreter', 'latex');
ylabel('$distance [m]$', 'interpreter', 'latex');
if save_plots
    saveas(gca, "E:/1. Unief/thesis/thesis_latex/figs/geom_diffs/distance.eps", 'epsc'); 
    saveas(gca, "E:/1. Unief/thesis/thesis_latex/figs/geom_diffs/distance.fig", 'fig'); 
end

