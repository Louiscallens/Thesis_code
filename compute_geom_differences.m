clear;

%% define variables
file_slack = "workspaces/thesis/final_problem3/high_accuracy_slack/iteration_";
iterations_slack = 1:10;
file_slack_post = "";

file_hp = "workspaces/thesis/final_problem3/high_accuracy_hp/iteration_";
iterations_hp = 1:5;
file_hp_post = "";

file_ref = "E:\1. Unief\thesis\thesis_latex\figs\final_problem3\reference_N_experiment\N_";
iterations_ref = 100:100:1500;
file_ref_post = "meta_data";

file2 = "E:\1. Unief\thesis\thesis_latex\figs\final_problem3\reference_N_experiment\N_1500meta_data";

%% perform computations
[diffs_slack, areas_slack, max_diffs_slack, svals_all_slack, fvals_all_slack] = get_geom_metrics(file_slack, iterations_slack, file_slack_post, file2);
[diffs_hp, areas_hp, max_diffs_hp, svals_all_hp, fvals_all_hp] = get_geom_metrics(file_hp, iterations_hp, file_hp_post, file2);
[diffs_ref, areas_ref, max_diffs_ref, svals_all_ref, fvals_all_ref] = get_geom_metrics(file_ref, iterations_ref, file_ref_post, file2);

save("workspaces/thesis/final_problem3/geom_diffs");

%% helper functions
function [diffs, areas, max_diffs, svals_all, fvals_all] = get_geom_metrics(file1, iterations, file1_post, file2)
    load(file2)
    res2 = res; M2 = M;
    
    diffs = [];
    areas = [];
    max_diffs = [];
    svals_all = {};
    fvals_all = {};
    
    for i = iterations
        file = file1 + num2str(i) + file1_post;
        load(file);
        res1 = res; M1 = M;

        %diffs = [diffs, get_geometric_diff(res1, M1, res2, M2, problem)];
        areas = [areas, get_geometric_area(res1, M1, res2, M2, problem)];
        [dist, svals, fvals] = get_max_euclidian_dist(res1, M1, res2, M2, problem);
        max_diffs = [max_diffs, dist];
        svals_all{i} = svals;
        fvals_all{i} = fvals;
    end
end

function diff = get_geometric_diff(res, M, res_ref, M_ref, problem)
    diff = integral(@(s) norm(evaluate_trajectory_e(res, M, s) - evaluate_trajectory_e(res_ref, M_ref, s)), 0, problem.myTrack.total_length, 'arrayValued', true);
end

function e = evaluate_trajectory_e(res, M, s)
    e = evaluate_trajectory(res, M, s);
    e = e(1);
end

function diff = get_geometric_area(res, M, res_ref, M_ref, problem)
    diff = integral(@(s) abs(evaluate_trajectory_e(res, M, s) - evaluate_trajectory_e(res_ref, M_ref, s)), 0, problem.myTrack.total_length, 'arrayValued', true);
end

function [dist, svals, fvals] = get_max_euclidian_dist(res1, M1, res2, M2, problem)
    svals = linspace(0, problem.myTrack.total_length, 10000);
    f = @(s) abs(evaluate_trajectory_e(res1, M1, s) - evaluate_trajectory_e(res2, M2, s));
    fvals = NaN + zeros(size(svals));
    for i = 1:length(svals)
        fvals(i) = f(svals(i));
    end
    dist = max(fvals);
end