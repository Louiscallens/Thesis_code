import casadi.*
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%% set up the problem
problem_switch = 0;
problem = setup_problem(problem_switch);

%% specify method parameters
method.maxIter = 1;
method.Nmin = 5;

method.save_plots = false;
method.plot_name = "";
method.og_plot_name = method.plot_name;

%% solve the problem

% initialize mesh
N = 15;
M = mesh(N, problem.myTrack.total_length, method.Nmin);
usedMeshes = {};

% initialize loop
converged = false;
iterCount = 1;
qualityMetrics = [];
res_previous = struct('X', {{NaN + zeros(problem.nx,0)}}, 'U', {{NaN + zeros(problem.nu,0)}});
M_previous = M;

% do the loop
while ~converged && iterCount <= method.maxIter
    % solve ocp
    res = solve_ocp(M, problem, M_previous, res_previous);
    
    % store some intermediate results
    usedMeshes{end+1} = M;
    [~, specifics] = get_quality_metric(res, M, problem.rhs, problem);
    quality_metrics(:,IterCount) = specifics;
    
    % do some plotting
    displayTrajectoryX(res, M, problem, method.save_plots, method.plot_name);
    displayTrajectoryU(res, M, problem, method.save_plots, method.plot_name);
    
    % check for convergence
    if iterCount == method.maxIter
        break;
    end
    
    % update the mesh
    M = get_new_mesh(res, M, problem, method);
    res_previous = res;
    iterCount = iterCount + 1;
end

%% display some final results
displayQualityMetrics(quality_metrics, method)