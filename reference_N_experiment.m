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
N_values = [100, 200, 300, 400, 500, 600];
exp_timings_total = [];
exp_timings_comp  = [];
exp_timings_build = [];
exp_iterCounts = [];
exp_tfs = [];
exp_nb_vars = [];

clear;
load('workspaces/thesis/final_problem3/reference_N_experiment/end_of_experiment');
nn_init = length(N_values);
N_values = [N_values, 1300, 1400, 1500];

for nn = nn_init+1:length(N_values)
    
method.method_select = 0; % 0: slackness-based method - 1: basic hp (patterson)
method.N = N_values(nn);
method.maxIter = 1;
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

method.save_plots = true;
method.save_every_iteration = false;
method.plot_iterations = [];
method.plot_name = "E:/1. Unief/thesis/thesis_latex/figs/final_problem3/reference_N_experiment/N_"+num2str(method.N);
method.og_plot_name = method.plot_name;
method.save_workspace = method.save_plots || method.save_every_iteration || ~isempty(method.plot_iterations);
method.workspace_name = "workspaces/thesis/final_problem3/reference_N_experiment/N_"+num2str(method.N);
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

t1 = tic;
[res, timing, iter_count, exit_code] = solve_ocp(M, problem, problem_switch, method, M_previous, res_previous);
a = toc(t1);
exp_timings_total(nn) = a;
exp_timings_comp(nn) = timing;
exp_timings_build(nn) = a - timing;
exp_iterCounts(nn) = iter_count;
exp_tfs(nn) = res.tf;
exp_nb_vars(nn) = count_nb_vars(M, problem);

save(method.og_plot_name+"meta_data");

end

save('workspaces/thesis/final_problem3/reference_N_experiment/end_of_experiment');