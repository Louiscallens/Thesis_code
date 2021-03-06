import casadi.*
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize','factory');

%% set up the problem
% 0: chicane - 1: smooth sine - 2: hairpin - 3: generic - 4: smooth hairpin - 5: circle - 6: zoomed chicane - 7: straight line
% 8: corner-cutting track - 9: smooth controls necessity - 10: masking - 11: spa
% effect
problem_switch = 11;
problem = setup_problem(problem_switch);

%% specify method parameters
method.method_select = 1; % 0: slackness-based method - 1: basic hp (patterson)
method.N = 40;
method.maxIter = 10;
method.Nmin = 5;
method.Nstep = 4;
method.Nmax = 12;
method.minUDegree = 0 + 2*method.method_select; %0: piecewise constant - 1: piecewise linear - 2: polynomial (control value for every collocation point)
if method.method_select == 0
    method.slack_performance_treshold = 1.0e-2;%1.0e-2;
    %method.slack_path_treshold = 0.5;%1.0e-1;%1.0e-2;
    method.slack_path_treshold = 1.0e-35;%1.0e-1;%1.0e-2;
    %method.err_treshold = 1.0e-1; %1.0e-4;
    method.err_treshold = 10^(-5); %1.0e-4;
    method.err_priority_treshold = 1.0;
    method.err_order_magnitude_treshold = 2;
else
    method.slack_performance_treshold = 1.0e30; method.slack_path_treshold = -1;
    method.err_treshold = 10^(-5); method.err_priority_treshold = 1.0e5;
    method.err_order_magnitude_treshold = inf;
end
method.minimal_interval_width = 1.0e-0;
method.use_viol_vars = true || method.method_select;
method.viol_cost_weight = method.use_viol_vars*1.0e10;
method.regularization_weight = 0*1.0e-2;%1.0e-4;

method.use_warm_start = true;

method.save_plots = false;
method.save_every_iteration = false;
method.plot_iterations = [];
method.plot_name = "E:/1. Unief/thesis/thesis_latex/figs/final_problem3/reference/N_600";
method.og_plot_name = method.plot_name;
method.save_workspace = method.save_plots || method.save_every_iteration || ~isempty(method.plot_iterations);
method.workspace_name = "workspaces/thesis/final_problem3/reference/N_600";
method.og_workspace_name = method.workspace_name;
method.skip_plot_position = method.save_plots || method.save_every_iteration || ~isempty(method.plot_iterations);
method.plot_metrics_separately = method.save_plots || method.save_every_iteration || ~isempty(method.plot_iterations);
if method.save_plots || method.save_every_iteration || ~isempty(method.plot_iterations)
    set(groot,'defaultAxesFontSize',14);
end

method.load_reference = false || (problem_switch == 6);
method.save_result = false;

%% solve the problem

% initialize mesh
M = mesh(method.N, problem.myTrack.total_length, method.Nmin, method.minUDegree, problem.nu, problem.disconts);
%M = maskingMesh2(method.N, problem.myTrack.total_length, method.Nmin, method.minUDegree, problem.nu, problem.disconts);
%if problem_switch == 10; M = maskingMesh(method.N, problem.myTrack.total_length, method.Nmin, method.minUDegree, problem.nu, problem.disconts); end
%load('mesh_chicane_trouble.mat');
%load('mesh_chicane_oscillations.mat');
usedMeshes = {};

% initialize loop
converged = false;
iterCount = 1;
qualityMetrics = [];
if method.load_reference
    %{
    if problem_switch == 6
        load(problem.reference_name_full);
        res_ref = res_previous;
        M_ref = M_previous;
        M_previous.s = M_previous.s(problem.N_first:problem.N_end);
        M_previous.Nk = M_previous.Nk(problem.N_first:problem.N_end);
        M_previous.Nu = M_previous.Nu(problem.N_first:problem.N_end);
        M_previous.sc = {};
        res = struct(); res.X = {}; res.U = {}; res.Yx = {}; res.Yu = {}; res.tc = {}; res.t = [];
        for k = problem.N_first:problem.N_end-1
            res.X{k+1-problem.N_first}  = res_ref.X{k};% + 1.0e-1.*rand(size(res_ref.X{k}));
            res.U{k+1-problem.N_first}  = res_ref.U{k};
            M_previous.sc{k+1-problem.N_first} = M_ref.sc{k};
        end
        res.X{problem.N_end+1-problem.N_first}  = res_ref.X{problem.N_end};
        res.U{problem.N_end+1-problem.N_first}  = res_ref.U{problem.N_end};
        res_previous = res;
    else
    %}
        load(problem.reference_name_full);
    %end
else
    res_previous = struct('X', {{NaN + zeros(problem.nx,0)}}, 'U', {{NaN + zeros(problem.nu,0)}});
    M_previous = M;
end
timings = [];
iterCounts = [];

%clear;
%load('workspaces/thesis/final_problem3/high_accuracy_hp/iteration_5');

% do the loop
while ~converged && iterCount <= method.maxIter
    % solve ocp
    [res, timing, iter_count, exit_code] = solve_ocp(M, problem, problem_switch, method, M_previous, res_previous);
    timings(iterCount) = timing;
    iterCounts(iterCount) = iter_count;
    if exit_code ~= 0
        break;
    end
    
    % store some intermediate results
    usedMeshes{end+1} = M;
    [~, specifics] = get_quality_metric(res, M, problem.rhs, problem, method);
    qualityMetrics(:,iterCount) = specifics;
    
    % do some plotting
    if method.save_every_iteration || ~isempty(method.plot_iterations)
        method.plot_name = method.og_plot_name + num2str(iterCount);
    end
    displayTrajectoryX(res, M, problem, method.save_plots && iterCount == method.maxIter || method.save_every_iteration || ismember(iterCount, method.plot_iterations), method.plot_name, method);
    displayTrajectoryU(res, M, problem, method.save_plots && iterCount == method.maxIter || method.save_every_iteration || ismember(iterCount, method.plot_iterations), method.plot_name, method);
    if method.plot_metrics_separately || problem_switch == 7 || problem_switch == 8 || problem_switch == 9
        displayQualityMetricsSeparately(qualityMetrics, method.save_plots && iterCount == method.maxIter || method.save_every_iteration || ismember(iterCount, method.plot_iterations), method.plot_name, method)
    else
        displayQualityMetrics(qualityMetrics, method.save_plots && iterCount == method.maxIter || method.save_every_iteration || ismember(iterCount, method.plot_iterations), method.plot_name, method)
    end
    [~, errs] = get_error_est(res, M, problem.rhs, method, problem.disconts);
    
    method.workspace_name = method.og_workspace_name + num2str(iterCount);
    save(method.workspace_name);
    
    % check for convergence
    if iterCount == method.maxIter
        break;
    end
    
    % update the mesh
    if method.use_warm_start
        M_previous = M;
        res_previous = res;
    end
    [M, updated] = get_new_mesh(res, M, errs, problem, method);
    if ~updated
        displayTrajectoryX(res, M, problem, method.save_plots, method.plot_name, method);
        displayTrajectoryU(res, M, problem, method.save_plots, method.plot_name, method);
        if method.plot_metrics_separately || problem_switch == 7 || problem_switch == 8 || problem_switch == 9
            displayQualityMetricsSeparately(qualityMetrics, method.save_plots, method.plot_name, method)
        else
            displayQualityMetrics(qualityMetrics, method.save_plots, method.plot_name, method)
        end
        break;
    end
    drawnow();
    iterCount = iterCount + 1;
    disp("-------- STARTING ITERATION "+num2str(iterCount)+" --------");
end

if method.save_result
    res_previous = res; M_previous = M;
    save(problem.reference_name+"_N_"+num2str(N)+".mat", 'res_previous', 'M_previous');
end

save(method.og_plot_name+"meta_data");

% list of figures:
% 1.  all states
% 2.  trajectory
% 3.  inputs
% 4.  quality metrics (end)
% 5.  relative errors
% 6.  mesh updates
% 7.  final time
% 8.  nb vars
% 9.  accuracy
% 10. feasibility