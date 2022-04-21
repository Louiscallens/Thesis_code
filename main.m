import casadi.*
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%% set up the problem
problem_switch = 6; % 0: chicane - 1: smooth sine - 2: hairpin - 3: generic - 4: smooth hairpin - 5: circle - 6: zoomed chicane
problem = setup_problem(problem_switch);

%% specify method parameters
method.N = 30;
method.maxIter = 2;
method.Nmin = 5;
method.Nstep = 5;
method.Nmax = 20;
method.minUDegree = 0; %0: piecewise constant - 1: piecewise linear - 2: polynomial (control value for every collocation point
method.slack_performance_treshold = 1.0e-2;%1.0e-2;
method.slack_path_treshold = 1.0e-10;%1.0e-2;
method.err_treshold = 1.0e-8;
method.err_priority_treshold = 1.0;
%method.slack_performance_treshold = 1.0e30; method.slack_path_treshold = -1;
%method.err_treshold = 1.0e-8; method.err_priority_treshold = 1.0e30;

method.save_plots = false;
method.plot_name = "figs/poster/hairpin";
method.og_plot_name = method.plot_name;

method.load_reference = true;
method.save_result = false;

%% solve the problem

% initialize mesh
M = mesh(method.N, problem.myTrack.total_length, method.Nmin, method.minUDegree, problem.disconts);
%load('mesh_chicane_trouble.mat');
%load('mesh_chicane_oscillations.mat');
usedMeshes = {};

% initialize loop
converged = false;
iterCount = 1;
qualityMetrics = [];
if method.load_reference
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
        load(problem.reference_name_full);
    end
else
    res_previous = struct('X', {{NaN + zeros(problem.nx,0)}}, 'U', {{NaN + zeros(problem.nu,0)}});
    M_previous = M;
end

% do the loop
while ~converged && iterCount <= method.maxIter
    % solve ocp
    res = solve_ocp(M, problem, M_previous, res_previous);
    
    % store some intermediate results
    usedMeshes{end+1} = M;
    [~, specifics] = get_quality_metric(res, M, problem.rhs, problem, method);
    qualityMetrics(:,iterCount) = specifics;
    
    % do some plotting
    displayTrajectoryX(res, M, problem, method.save_plots && iterCount == method.maxIter, method.plot_name);
    displayTrajectoryU(res, M, problem, method.save_plots && iterCount == method.maxIter, method.plot_name);
    displayQualityMetrics(qualityMetrics, method.save_plots && iterCount == method.maxIter, method.plot_name)
    
    % check for convergence
    if iterCount == method.maxIter
        break;
    end
    
    % update the mesh
    %M_previous = M;
    %res_previous = res;
    M = get_new_mesh(res, M, problem, method);
    iterCount = iterCount + 1;
    disp("-------- STARTING ITERATION "+num2str(iterCount)+" --------");
end

if method.save_result
    res_previous = res; M_previous = M;
    save(problem.reference_name+"_N_"+num2str(N)+".mat", 'res_previous', 'M_previous');
end


% list of figures:
% 1. all states
% 2. trajectory
% 3. inputs
% 4. quality metrics (end)
% 5. relative errors
% 6. mesh updates