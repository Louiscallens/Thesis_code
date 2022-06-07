load('workspaces/thesis/final_problem3/reference_N_experiment/end_of_experiment');
set(groot,'defaultAxesFontSize',14);

save_figs = false;

figure;
plot(N_values, exp_timings_total, '.-k', 'linewidth', 1); hold on;
plot(N_values, exp_timings_comp,  '.-b', 'linewidth', 1);
plot(N_values, exp_timings_build, '.-r', 'linewidth', 1);
xlim([N_values(1), N_values(end)]);
xlabel('$N$', 'interpreter', 'latex');
ylabel('$time [s]$', 'interpreter', 'latex');
grid on;
if save_figs; saveas(gca, "E:/1. Unief/thesis/thesis_latex/figs/final_problem3/reference_N_experiment/timings.eps", 'epsc'); end

figure;
plot(N_values, exp_iterCounts, '.-b', 'linewidth', 1);
xlabel('$N$', 'interpreter', 'latex');
ylabel('number of iterations', 'interpreter', 'latex');
xlim([N_values(1), N_values(end)]);
ylim([40, 60]);
grid on;
if save_figs; saveas(gca, "E:/1. Unief/thesis/thesis_latex/figs/final_problem3/reference_N_experiment/iterCounts.eps", 'epsc'); end

figure;
plot(N_values, exp_tfs, '.-b', 'linewidth', 1);
xlabel('$N$', 'interpreter', 'latex');
ylabel('$t_f [s]$', 'interpreter', 'latex');
xlim([N_values(1), N_values(end)]);
grid on;
if save_figs; saveas(gca, "E:/1. Unief/thesis/thesis_latex/figs/final_problem3/reference_N_experiment/tfs.eps", 'epsc'); end

figure;
plot(N_values, exp_nb_vars, '.-b', 'linewidth', 1);
xlabel('$N$', 'interpreter', 'latex');
ylabel('number of variables', 'interpreter', 'latex');
xlim([N_values(1), N_values(end)]);
grid on;
if save_figs; saveas(gca, "E:/1. Unief/thesis/thesis_latex/figs/final_problem3/reference_N_experiment/nbVars.eps", 'epsc'); end
